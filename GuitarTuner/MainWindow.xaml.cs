using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Reflection.Emit;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using System.Timers;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Forms;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using FftSharp;
using FftSharp.Windows;
using MathNet.Numerics;
using NAudio.Wave;
using ScottPlot.Plottable;

namespace GuitarTuner
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : System.Windows.Window
    {
        //const int SAMPLE_RATE = 44_100;
        const int SAMPLE_RATE = 48_000;

        const int BIT_DEPTH = 16;

        const int CHANNEL_COUNT = 1;

        const int BUFFER_MS = 420;

        const int NUM_HPS = 5;

        //const int BUFFER_MS = 110;

        VLine peakLine;

        /// <summary>
        /// The audio values to be plotted
        /// </summary>
        readonly double[] AudioValues;

        readonly double[] FftValues;

        System.Windows.Forms.Timer timer = new System.Windows.Forms.Timer();

        public MainWindow()
        {
            InitializeComponent();

            AudioValues = new double[SAMPLE_RATE * BUFFER_MS / 1000];
            double[] paddedAudio = Pad.ZeroPad(AudioValues);
            double[] fftMag = FftSharp.Transform.FFTmagnitude(paddedAudio);
            FftValues = new double[fftMag.Length];

            double fftPeriod = FftSharp.Transform.FFTfreqPeriod(SAMPLE_RATE, fftMag.Length);

            PlotGraph.Plot.AddSignal(AudioValues, SAMPLE_RATE / 1000);
            PlotGraph.Plot.YLabel("Level");
            PlotGraph.Plot.XLabel("Time (milliseconds)");
            PlotGraph.Refresh();

            var rate = 2.0 * fftMag.Length / SAMPLE_RATE * NUM_HPS;
            FrequencyGraph.Plot.AddSignal(FftValues, rate);
            FrequencyGraph.Plot.YLabel("Spectral Power");
            FrequencyGraph.Plot.XLabel("Frequency (kHz)");
            FrequencyGraph.Refresh();

            timer.Interval = BUFFER_MS / 2;
            timer.Tick += Timer_Tick;
            timer.Tick += Timer_Tick_Frequency;
            timer.Start();

            BeginRecording();
        }

        private void BeginRecording()
        {
            var waveIn = new WaveInEvent
            {
                DeviceNumber = 0, // indicates which microphone to use
                WaveFormat = new WaveFormat(SAMPLE_RATE, BIT_DEPTH, CHANNEL_COUNT),
                BufferMilliseconds = BUFFER_MS
            };
            waveIn.DataAvailable += WaveIn_DataAvailable2;
            waveIn.StartRecording();
        }

        private void Timer_Tick(object? sender, EventArgs e)
        {
            int level = (int)AudioValues.Max();

            // auto-scale the maximum progressbar level
            if (level > ProgressBar.Maximum)
            {
                ProgressBar.Maximum = level;
            }
            ProgressBar.Value = level;

            // auto-scale the plot Y axis limits
            var currentLimits = PlotGraph.Plot.GetAxisLimits();
            PlotGraph.Plot.SetAxisLimits(
                yMin: Math.Min(currentLimits.YMin, -level),
                yMax: Math.Max(currentLimits.YMax, level)
            );

            // request a redraw using a non-blocking render queue
            PlotGraph.RefreshRequest();
        }

        private void Timer_Tick_Frequency(object? sender, EventArgs e)
        {
            var window = new Hanning();
            double[] windowed = window.Apply(AudioValues);
            double[] paddedAudio = Pad.ZeroPad(windowed);
            double[] filtered = Filter.HighPass(paddedAudio, SAMPLE_RATE, 62);
            double[] fftMag = FftSharp.Transform.FFTmagnitude(filtered);
            fftMag = HarmonicProductSpectrum(InterpolateAndNormalize(fftMag));
            double[] fftFreq = FftSharp.Transform.FFTfreq(SAMPLE_RATE, fftMag.Length);

            // find the frequency peak
            int peakIndex = 0;
            for (int i = 0; i < fftMag.Length; i++)
            {
                if (fftMag[i] > fftMag[peakIndex])
                {
                    peakIndex = i;
                }
            }
            double peakPower = fftMag[peakIndex];
            double peakFrequency = fftFreq[peakIndex] / 5.0;
            Array.Copy(fftMag, FftValues, fftMag.Length);
            var closestNote = FindClosestNote(peakFrequency);
            DetectedFrequencyLabel.Content = $"Fundamental Frequency: {peakFrequency:N2} Hz";
            NearestFrequencyLabel.Content = $"{closestNote.note} {closestNote.pitch:N2}Hz";

            FrequencyGraph.Plot.SetAxisLimits(
                xMin: 0,
                xMax: 1000,
                yMin: 0,
                yMax: peakPower + (peakPower * 0.1)
            //yMax: Math.Max(peakPower, FrequencyGraph.Plot.GetAxisLimits().YMax)
            );

            // request a redraw using a non-blocking render queue
            FrequencyGraph.RefreshRequest();
        }

        private double[] InterpolateAndNormalize(double[] fftMag)
        {
            // Interpolate
            var needToInterpolate = new List<double>();
            var xPoints = Enumerable.Range(0, fftMag.Length).Select(x => (double)x).ToArray();
            double lastInserted = 0;
            while (lastInserted < (fftMag.Length - 0.01))
            {
                needToInterpolate.Add(lastInserted);
                lastInserted += 1 / 5.0;
            }
            var interpolation = Interpolate.Linear(xPoints, fftMag);
            var interpolatedMagnitudes = needToInterpolate.Select(interpolation.Interpolate);

            // normalize
            var ratio = 1.0 / interpolatedMagnitudes.Max();
            return interpolatedMagnitudes.Select(x => x * ratio).ToArray();
        }

        public double[] HarmonicProductSpectrum(double[] data)
        {
            double[] hps2 = Downsample(data, 2);
            double[] hps3 = Downsample(data, 3);
            double[] hps4 = Downsample(data, 4);
            double[] hps5 = Downsample(data, 5);
            double[] array = new double[hps5.Length];

            for (int i = 0; i < array.Length; i++)
            {
                checked
                {
                    array[i] = data[i] * hps2[i] * hps3[i] * hps4[i] * hps5[i];
                }
            }
            return array;
        }

        public double[] Downsample(double[] data, int n)
        {
            double[] array = new double[Convert.ToInt32(Math.Ceiling(data.Length * 1.0 / n))];
            for (int i = 0; i < array.Length; i++)
            {
                array[i] = data[i * n];
            }
            return array;
        }

        List<string> allNotes = new List<string> { "A", "A#", "B", "C", "C#", "D", "D#", "E", "F", "F#", "G", "G#" };
        const double CONCERT_PITCH = 440d;
        private (string note, double pitch) FindClosestNote(double pitch)
        {
            var rawIndex = (int)Math.Round(Math.Log2(pitch / CONCERT_PITCH) * 12);
            var adjustedIndex = rawIndex;
            while (adjustedIndex < 0)
            {
                adjustedIndex += 12;
            }
            var closestNote = allNotes[adjustedIndex % 12] + (4 + Math.Floor((rawIndex + 9) / 12.0)).ToString();
            var closestPitch = CONCERT_PITCH * Math.Pow(2, rawIndex / 12.0);
            return (closestNote, closestPitch);
        }

        public void WaveIn_DataAvailable2(object? sender, WaveInEventArgs e)
        {
            for (int i = 0; i < e.Buffer.Length / 2; i++)
            {
                AudioValues[i] = BitConverter.ToInt16(e.Buffer, i * 2);
            }
        }

    }
}
