using System;
using System.Collections.Generic;
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
using NAudio.Wave;

namespace GuitarTuner
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : System.Windows.Window
    {
        const int SAMPLE_RATE = 44_100;

        const int BIT_DEPTH = 16;

        const int CHANNEL_COUNT = 1;

        const int BUFFER_MS = 20;

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

            FrequencyGraph.Plot.AddSignal(FftValues, fftPeriod);
            FrequencyGraph.Plot.YLabel("Spectral Power");
            FrequencyGraph.Plot.XLabel("Frequency (kHz)");
            FrequencyGraph.Refresh();

            timer.Interval = BUFFER_MS;
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
            double[] paddedAudio = FftSharp.Pad.ZeroPad(AudioValues);
            double[] fftMag = FftSharp.Transform.FFTpower(paddedAudio);
            Array.Copy(fftMag, FftValues, fftMag.Length);

            // find the frequency peak
            int peakIndex = 0;
            for (int i = 0; i < fftMag.Length; i++)
            {
                if (fftMag[i] > fftMag[peakIndex])
                    peakIndex = i;
            }
            double fftPeriod = FftSharp.Transform.FFTfreqPeriod(SAMPLE_RATE, fftMag.Length);
            double peakFrequency = fftPeriod * peakIndex;
            //label1.Text = $"Peak Frequency: {peakFrequency:N0} Hz";

            // auto-scale the plot Y axis limits
            double fftPeakMag = fftMag.Max();
            double plotYMax = FrequencyGraph.Plot.GetAxisLimits().YMax;
            FrequencyGraph.Plot.SetAxisLimits(
                xMin: 0,
                xMax: 3,
                yMin: 0,
                yMax: Math.Max(fftPeakMag, plotYMax));

            // request a redraw using a non-blocking render queue
            FrequencyGraph.RefreshRequest();
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
