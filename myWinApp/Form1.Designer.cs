namespace myWinApp
{
    partial class Form1
    {
        /// <summary>
        /// 設計工具所需的變數。
        /// </summary>
        private System.ComponentModel.IContainer components = null;

        /// <summary>
        /// 清除任何使用中的資源。
        /// </summary>
        /// <param name="disposing">如果應該處置 Managed 資源則為 true，否則為 false。</param>
        protected override void Dispose(bool disposing)
        {
            if (disposing && (components != null))
            {
                components.Dispose();
            }
            base.Dispose(disposing);
        }

        #region Windows Form 設計工具產生的程式碼

        /// <summary>
        /// 此為設計工具支援所需的方法 - 請勿使用程式碼編輯器修改
        /// 這個方法的內容。
        /// </summary>
        private void InitializeComponent()
        {
            System.Windows.Forms.DataVisualization.Charting.ChartArea chartArea4 = new System.Windows.Forms.DataVisualization.Charting.ChartArea();
            System.Windows.Forms.DataVisualization.Charting.Series series7 = new System.Windows.Forms.DataVisualization.Charting.Series();
            System.Windows.Forms.DataVisualization.Charting.Series series8 = new System.Windows.Forms.DataVisualization.Charting.Series();
            this.Simulation = new System.Windows.Forms.TabControl();
            this.priceTest2 = new System.Windows.Forms.TabPage();
            this.pricingResult = new System.Windows.Forms.GroupBox();
            this.textBox_putPrice = new System.Windows.Forms.TextBox();
            this.textBox_callPrice = new System.Windows.Forms.TextBox();
            this.label12 = new System.Windows.Forms.Label();
            this.label13 = new System.Windows.Forms.Label();
            this.simulationPara = new System.Windows.Forms.GroupBox();
            this.textBox_seed = new System.Windows.Forms.TextBox();
            this.label15 = new System.Windows.Forms.Label();
            this.label14 = new System.Windows.Forms.Label();
            this.textBox_pathCnt = new System.Windows.Forms.TextBox();
            this.button2 = new System.Windows.Forms.Button();
            this.msgGroupBox = new System.Windows.Forms.GroupBox();
            this.clear = new System.Windows.Forms.Button();
            this.msgBox = new System.Windows.Forms.RichTextBox();
            this.hestonPara = new System.Windows.Forms.GroupBox();
            this.textBox_sigma = new System.Windows.Forms.TextBox();
            this.label7 = new System.Windows.Forms.Label();
            this.label8 = new System.Windows.Forms.Label();
            this.textBox_theta = new System.Windows.Forms.TextBox();
            this.label9 = new System.Windows.Forms.Label();
            this.label10 = new System.Windows.Forms.Label();
            this.textBox_kappa = new System.Windows.Forms.TextBox();
            this.label11 = new System.Windows.Forms.Label();
            this.textBox_rho = new System.Windows.Forms.TextBox();
            this.optionPara = new System.Windows.Forms.GroupBox();
            this.textBox_rf = new System.Windows.Forms.TextBox();
            this.textBox_T = new System.Windows.Forms.TextBox();
            this.label6 = new System.Windows.Forms.Label();
            this.label5 = new System.Windows.Forms.Label();
            this.textBox_var0 = new System.Windows.Forms.TextBox();
            this.label4 = new System.Windows.Forms.Label();
            this.label3 = new System.Windows.Forms.Label();
            this.textBox_k = new System.Windows.Forms.TextBox();
            this.label2 = new System.Windows.Forms.Label();
            this.textBox_s0 = new System.Windows.Forms.TextBox();
            this.tabPage1 = new System.Windows.Forms.TabPage();
            this.textBox_rho_old = new System.Windows.Forms.TextBox();
            this.label1 = new System.Windows.Forms.Label();
            this.richTextBox1 = new System.Windows.Forms.RichTextBox();
            this.button1 = new System.Windows.Forms.Button();
            this.chart1 = new System.Windows.Forms.DataVisualization.Charting.Chart();
            this.priceTest1 = new System.Windows.Forms.TabPage();
            this.priceResult = new System.Windows.Forms.TextBox();
            this.listBox1 = new System.Windows.Forms.ListBox();
            this.pricingButton = new System.Windows.Forms.Button();
            this.Simulation.SuspendLayout();
            this.priceTest2.SuspendLayout();
            this.pricingResult.SuspendLayout();
            this.simulationPara.SuspendLayout();
            this.msgGroupBox.SuspendLayout();
            this.hestonPara.SuspendLayout();
            this.optionPara.SuspendLayout();
            this.tabPage1.SuspendLayout();
            ((System.ComponentModel.ISupportInitialize)(this.chart1)).BeginInit();
            this.priceTest1.SuspendLayout();
            this.SuspendLayout();
            // 
            // Simulation
            // 
            this.Simulation.Controls.Add(this.priceTest2);
            this.Simulation.Controls.Add(this.tabPage1);
            this.Simulation.Controls.Add(this.priceTest1);
            this.Simulation.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.Simulation.Location = new System.Drawing.Point(0, 0);
            this.Simulation.Name = "Simulation";
            this.Simulation.SelectedIndex = 0;
            this.Simulation.Size = new System.Drawing.Size(1800, 798);
            this.Simulation.TabIndex = 6;
            // 
            // priceTest2
            // 
            this.priceTest2.BackColor = System.Drawing.Color.AliceBlue;
            this.priceTest2.Controls.Add(this.pricingResult);
            this.priceTest2.Controls.Add(this.simulationPara);
            this.priceTest2.Controls.Add(this.msgGroupBox);
            this.priceTest2.Controls.Add(this.hestonPara);
            this.priceTest2.Controls.Add(this.optionPara);
            this.priceTest2.Location = new System.Drawing.Point(8, 47);
            this.priceTest2.Name = "priceTest2";
            this.priceTest2.Size = new System.Drawing.Size(1784, 743);
            this.priceTest2.TabIndex = 2;
            this.priceTest2.Text = "priceTest2";
            // 
            // pricingResult
            // 
            this.pricingResult.Controls.Add(this.textBox_putPrice);
            this.pricingResult.Controls.Add(this.textBox_callPrice);
            this.pricingResult.Controls.Add(this.label12);
            this.pricingResult.Controls.Add(this.label13);
            this.pricingResult.Location = new System.Drawing.Point(930, 313);
            this.pricingResult.Name = "pricingResult";
            this.pricingResult.Size = new System.Drawing.Size(392, 182);
            this.pricingResult.TabIndex = 13;
            this.pricingResult.TabStop = false;
            this.pricingResult.Text = "pricing result";
            // 
            // textBox_putPrice
            // 
            this.textBox_putPrice.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.textBox_putPrice.Location = new System.Drawing.Point(150, 119);
            this.textBox_putPrice.Name = "textBox_putPrice";
            this.textBox_putPrice.Size = new System.Drawing.Size(226, 46);
            this.textBox_putPrice.TabIndex = 10;
            this.textBox_putPrice.Text = "0";
            // 
            // textBox_callPrice
            // 
            this.textBox_callPrice.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.textBox_callPrice.Location = new System.Drawing.Point(150, 54);
            this.textBox_callPrice.Name = "textBox_callPrice";
            this.textBox_callPrice.Size = new System.Drawing.Size(226, 46);
            this.textBox_callPrice.TabIndex = 9;
            this.textBox_callPrice.Text = "0";
            // 
            // label12
            // 
            this.label12.AutoSize = true;
            this.label12.Location = new System.Drawing.Point(32, 54);
            this.label12.Name = "label12";
            this.label12.Size = new System.Drawing.Size(57, 32);
            this.label12.TabIndex = 9;
            this.label12.Text = "call";
            // 
            // label13
            // 
            this.label13.AutoSize = true;
            this.label13.Location = new System.Drawing.Point(32, 119);
            this.label13.Name = "label13";
            this.label13.Size = new System.Drawing.Size(53, 32);
            this.label13.TabIndex = 9;
            this.label13.Text = "put";
            // 
            // simulationPara
            // 
            this.simulationPara.Controls.Add(this.textBox_seed);
            this.simulationPara.Controls.Add(this.label15);
            this.simulationPara.Controls.Add(this.label14);
            this.simulationPara.Controls.Add(this.textBox_pathCnt);
            this.simulationPara.Controls.Add(this.button2);
            this.simulationPara.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.simulationPara.Location = new System.Drawing.Point(930, 33);
            this.simulationPara.Name = "simulationPara";
            this.simulationPara.Size = new System.Drawing.Size(392, 272);
            this.simulationPara.TabIndex = 12;
            this.simulationPara.TabStop = false;
            this.simulationPara.Text = "simulation parameter";
            // 
            // textBox_seed
            // 
            this.textBox_seed.Location = new System.Drawing.Point(150, 132);
            this.textBox_seed.Name = "textBox_seed";
            this.textBox_seed.Size = new System.Drawing.Size(226, 46);
            this.textBox_seed.TabIndex = 11;
            this.textBox_seed.Text = "1234";
            // 
            // label15
            // 
            this.label15.AutoSize = true;
            this.label15.Location = new System.Drawing.Point(32, 135);
            this.label15.Name = "label15";
            this.label15.Size = new System.Drawing.Size(68, 32);
            this.label15.TabIndex = 10;
            this.label15.Text = "seed";
            // 
            // label14
            // 
            this.label14.AutoSize = true;
            this.label14.Location = new System.Drawing.Point(6, 66);
            this.label14.Name = "label14";
            this.label14.Size = new System.Drawing.Size(122, 32);
            this.label14.TabIndex = 9;
            this.label14.Text = "# of path";
            // 
            // textBox_pathCnt
            // 
            this.textBox_pathCnt.Location = new System.Drawing.Point(150, 63);
            this.textBox_pathCnt.Name = "textBox_pathCnt";
            this.textBox_pathCnt.Size = new System.Drawing.Size(226, 46);
            this.textBox_pathCnt.TabIndex = 9;
            this.textBox_pathCnt.Text = "10000";
            // 
            // button2
            // 
            this.button2.Location = new System.Drawing.Point(12, 197);
            this.button2.Name = "button2";
            this.button2.Size = new System.Drawing.Size(364, 53);
            this.button2.TabIndex = 1;
            this.button2.Text = "price";
            this.button2.UseVisualStyleBackColor = true;
            this.button2.Click += new System.EventHandler(this.button2_Click);
            // 
            // msgGroupBox
            // 
            this.msgGroupBox.Controls.Add(this.clear);
            this.msgGroupBox.Controls.Add(this.msgBox);
            this.msgGroupBox.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.msgGroupBox.Location = new System.Drawing.Point(44, 501);
            this.msgGroupBox.Name = "msgGroupBox";
            this.msgGroupBox.Size = new System.Drawing.Size(1262, 205);
            this.msgGroupBox.TabIndex = 11;
            this.msgGroupBox.TabStop = false;
            this.msgGroupBox.Text = "message";
            // 
            // clear
            // 
            this.clear.Location = new System.Drawing.Point(1122, 30);
            this.clear.Name = "clear";
            this.clear.Size = new System.Drawing.Size(134, 52);
            this.clear.TabIndex = 1;
            this.clear.Text = "clear";
            this.clear.UseVisualStyleBackColor = true;
            this.clear.Click += new System.EventHandler(this.clear_Click);
            // 
            // msgBox
            // 
            this.msgBox.Location = new System.Drawing.Point(29, 53);
            this.msgBox.Name = "msgBox";
            this.msgBox.Size = new System.Drawing.Size(1204, 138);
            this.msgBox.TabIndex = 0;
            this.msgBox.Text = "";
            // 
            // hestonPara
            // 
            this.hestonPara.Controls.Add(this.textBox_sigma);
            this.hestonPara.Controls.Add(this.label7);
            this.hestonPara.Controls.Add(this.label8);
            this.hestonPara.Controls.Add(this.textBox_theta);
            this.hestonPara.Controls.Add(this.label9);
            this.hestonPara.Controls.Add(this.label10);
            this.hestonPara.Controls.Add(this.textBox_kappa);
            this.hestonPara.Controls.Add(this.label11);
            this.hestonPara.Controls.Add(this.textBox_rho);
            this.hestonPara.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.hestonPara.Location = new System.Drawing.Point(486, 33);
            this.hestonPara.Name = "hestonPara";
            this.hestonPara.Size = new System.Drawing.Size(410, 462);
            this.hestonPara.TabIndex = 10;
            this.hestonPara.TabStop = false;
            this.hestonPara.Text = "heston parameters";
            // 
            // textBox_sigma
            // 
            this.textBox_sigma.Location = new System.Drawing.Point(167, 277);
            this.textBox_sigma.Name = "textBox_sigma";
            this.textBox_sigma.Size = new System.Drawing.Size(226, 46);
            this.textBox_sigma.TabIndex = 8;
            this.textBox_sigma.Text = "0.3";
            // 
            // label7
            // 
            this.label7.AutoSize = true;
            this.label7.Location = new System.Drawing.Point(38, 383);
            this.label7.Name = "label7";
            this.label7.Size = new System.Drawing.Size(0, 32);
            this.label7.TabIndex = 7;
            // 
            // label8
            // 
            this.label8.AutoSize = true;
            this.label8.Location = new System.Drawing.Point(38, 280);
            this.label8.Name = "label8";
            this.label8.Size = new System.Drawing.Size(86, 32);
            this.label8.TabIndex = 6;
            this.label8.Text = "sigma";
            // 
            // textBox_theta
            // 
            this.textBox_theta.Location = new System.Drawing.Point(167, 204);
            this.textBox_theta.Name = "textBox_theta";
            this.textBox_theta.Size = new System.Drawing.Size(226, 46);
            this.textBox_theta.TabIndex = 5;
            this.textBox_theta.Text = "0.04";
            // 
            // label9
            // 
            this.label9.AutoSize = true;
            this.label9.Location = new System.Drawing.Point(38, 207);
            this.label9.Name = "label9";
            this.label9.Size = new System.Drawing.Size(72, 32);
            this.label9.TabIndex = 4;
            this.label9.Text = "theta";
            // 
            // label10
            // 
            this.label10.AutoSize = true;
            this.label10.Location = new System.Drawing.Point(38, 135);
            this.label10.Name = "label10";
            this.label10.Size = new System.Drawing.Size(86, 32);
            this.label10.TabIndex = 3;
            this.label10.Text = "kappa";
            // 
            // textBox_kappa
            // 
            this.textBox_kappa.Location = new System.Drawing.Point(167, 132);
            this.textBox_kappa.Name = "textBox_kappa";
            this.textBox_kappa.Size = new System.Drawing.Size(226, 46);
            this.textBox_kappa.TabIndex = 2;
            this.textBox_kappa.Text = "1.5";
            // 
            // label11
            // 
            this.label11.AutoSize = true;
            this.label11.Location = new System.Drawing.Point(38, 63);
            this.label11.Name = "label11";
            this.label11.Size = new System.Drawing.Size(55, 32);
            this.label11.TabIndex = 1;
            this.label11.Text = "rho";
            // 
            // textBox_rho
            // 
            this.textBox_rho.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.textBox_rho.Location = new System.Drawing.Point(167, 60);
            this.textBox_rho.Name = "textBox_rho";
            this.textBox_rho.Size = new System.Drawing.Size(226, 46);
            this.textBox_rho.TabIndex = 0;
            this.textBox_rho.Text = "-0.9";
            // 
            // optionPara
            // 
            this.optionPara.Controls.Add(this.textBox_rf);
            this.optionPara.Controls.Add(this.textBox_T);
            this.optionPara.Controls.Add(this.label6);
            this.optionPara.Controls.Add(this.label5);
            this.optionPara.Controls.Add(this.textBox_var0);
            this.optionPara.Controls.Add(this.label4);
            this.optionPara.Controls.Add(this.label3);
            this.optionPara.Controls.Add(this.textBox_k);
            this.optionPara.Controls.Add(this.label2);
            this.optionPara.Controls.Add(this.textBox_s0);
            this.optionPara.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.optionPara.Location = new System.Drawing.Point(44, 33);
            this.optionPara.Name = "optionPara";
            this.optionPara.Size = new System.Drawing.Size(410, 462);
            this.optionPara.TabIndex = 0;
            this.optionPara.TabStop = false;
            this.optionPara.Text = "option parameters";
            // 
            // textBox_rf
            // 
            this.textBox_rf.Location = new System.Drawing.Point(167, 348);
            this.textBox_rf.Name = "textBox_rf";
            this.textBox_rf.Size = new System.Drawing.Size(226, 46);
            this.textBox_rf.TabIndex = 9;
            this.textBox_rf.Text = "0.001521";
            // 
            // textBox_T
            // 
            this.textBox_T.Location = new System.Drawing.Point(167, 277);
            this.textBox_T.Name = "textBox_T";
            this.textBox_T.Size = new System.Drawing.Size(226, 46);
            this.textBox_T.TabIndex = 8;
            this.textBox_T.Text = "1.0";
            // 
            // label6
            // 
            this.label6.AutoSize = true;
            this.label6.Location = new System.Drawing.Point(38, 351);
            this.label6.Name = "label6";
            this.label6.Size = new System.Drawing.Size(35, 32);
            this.label6.TabIndex = 7;
            this.label6.Text = "rf";
            // 
            // label5
            // 
            this.label5.AutoSize = true;
            this.label5.Location = new System.Drawing.Point(38, 280);
            this.label5.Name = "label5";
            this.label5.Size = new System.Drawing.Size(33, 32);
            this.label5.TabIndex = 6;
            this.label5.Text = "T";
            // 
            // textBox_var0
            // 
            this.textBox_var0.Location = new System.Drawing.Point(167, 204);
            this.textBox_var0.Name = "textBox_var0";
            this.textBox_var0.Size = new System.Drawing.Size(226, 46);
            this.textBox_var0.TabIndex = 5;
            this.textBox_var0.Text = "0.00770547621786487";
            // 
            // label4
            // 
            this.label4.AutoSize = true;
            this.label4.Location = new System.Drawing.Point(38, 207);
            this.label4.Name = "label4";
            this.label4.Size = new System.Drawing.Size(68, 32);
            this.label4.TabIndex = 4;
            this.label4.Text = "var0";
            // 
            // label3
            // 
            this.label3.AutoSize = true;
            this.label3.Location = new System.Drawing.Point(38, 135);
            this.label3.Name = "label3";
            this.label3.Size = new System.Drawing.Size(37, 32);
            this.label3.TabIndex = 3;
            this.label3.Text = "K";
            // 
            // textBox_k
            // 
            this.textBox_k.Location = new System.Drawing.Point(167, 132);
            this.textBox_k.Name = "textBox_k";
            this.textBox_k.Size = new System.Drawing.Size(226, 46);
            this.textBox_k.TabIndex = 2;
            this.textBox_k.Text = "100.0";
            // 
            // label2
            // 
            this.label2.AutoSize = true;
            this.label2.Location = new System.Drawing.Point(38, 63);
            this.label2.Name = "label2";
            this.label2.Size = new System.Drawing.Size(47, 32);
            this.label2.TabIndex = 1;
            this.label2.Text = "S0";
            // 
            // textBox_s0
            // 
            this.textBox_s0.Font = new System.Drawing.Font("新細明體", 12F, System.Drawing.FontStyle.Regular, System.Drawing.GraphicsUnit.Point, ((byte)(136)));
            this.textBox_s0.Location = new System.Drawing.Point(167, 60);
            this.textBox_s0.Name = "textBox_s0";
            this.textBox_s0.Size = new System.Drawing.Size(226, 46);
            this.textBox_s0.TabIndex = 0;
            this.textBox_s0.Text = "101.52";
            // 
            // tabPage1
            // 
            this.tabPage1.BackColor = System.Drawing.Color.LightCyan;
            this.tabPage1.BackgroundImageLayout = System.Windows.Forms.ImageLayout.None;
            this.tabPage1.Controls.Add(this.textBox_rho_old);
            this.tabPage1.Controls.Add(this.label1);
            this.tabPage1.Controls.Add(this.richTextBox1);
            this.tabPage1.Controls.Add(this.button1);
            this.tabPage1.Controls.Add(this.chart1);
            this.tabPage1.Location = new System.Drawing.Point(8, 47);
            this.tabPage1.Name = "tabPage1";
            this.tabPage1.Padding = new System.Windows.Forms.Padding(3);
            this.tabPage1.Size = new System.Drawing.Size(1784, 743);
            this.tabPage1.TabIndex = 0;
            this.tabPage1.Text = "volatility path";
            // 
            // textBox_rho_old
            // 
            this.textBox_rho_old.Location = new System.Drawing.Point(134, 127);
            this.textBox_rho_old.Name = "textBox_rho_old";
            this.textBox_rho_old.Size = new System.Drawing.Size(238, 46);
            this.textBox_rho_old.TabIndex = 9;
            // 
            // label1
            // 
            this.label1.AutoSize = true;
            this.label1.Font = new System.Drawing.Font("新細明體", 12F);
            this.label1.Location = new System.Drawing.Point(41, 131);
            this.label1.Name = "label1";
            this.label1.Size = new System.Drawing.Size(55, 32);
            this.label1.TabIndex = 10;
            this.label1.Text = "rho";
            // 
            // richTextBox1
            // 
            this.richTextBox1.Location = new System.Drawing.Point(134, 200);
            this.richTextBox1.Name = "richTextBox1";
            this.richTextBox1.Size = new System.Drawing.Size(238, 206);
            this.richTextBox1.TabIndex = 8;
            this.richTextBox1.Text = "";
            // 
            // button1
            // 
            this.button1.Location = new System.Drawing.Point(35, 200);
            this.button1.Name = "button1";
            this.button1.Size = new System.Drawing.Size(77, 48);
            this.button1.TabIndex = 7;
            this.button1.Text = "抽";
            this.button1.UseVisualStyleBackColor = true;
            // 
            // chart1
            // 
            this.chart1.BackColor = System.Drawing.Color.Transparent;
            chartArea4.AxisX.MajorGrid.Enabled = false;
            chartArea4.AxisX2.MajorGrid.Enabled = false;
            chartArea4.AxisY.Crossing = 1.7976931348623157E+308D;
            chartArea4.AxisY.LabelStyle.Format = "0.00";
            chartArea4.AxisY.LineColor = System.Drawing.Color.Transparent;
            chartArea4.AxisY.MajorGrid.Enabled = false;
            chartArea4.AxisY2.LabelStyle.Format = "0.00";
            chartArea4.AxisY2.MajorGrid.Enabled = false;
            chartArea4.AxisY2.TitleForeColor = System.Drawing.Color.White;
            chartArea4.BackColor = System.Drawing.Color.Transparent;
            chartArea4.Name = "ChartArea1";
            this.chart1.ChartAreas.Add(chartArea4);
            this.chart1.Location = new System.Drawing.Point(430, 123);
            this.chart1.Name = "chart1";
            series7.ChartArea = "ChartArea1";
            series7.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series7.Name = "volatility";
            series7.YAxisType = System.Windows.Forms.DataVisualization.Charting.AxisType.Secondary;
            series8.ChartArea = "ChartArea1";
            series8.ChartType = System.Windows.Forms.DataVisualization.Charting.SeriesChartType.Line;
            series8.Name = "stock";
            this.chart1.Series.Add(series7);
            this.chart1.Series.Add(series8);
            this.chart1.Size = new System.Drawing.Size(837, 401);
            this.chart1.TabIndex = 6;
            this.chart1.Text = "B";
            // 
            // priceTest1
            // 
            this.priceTest1.BackColor = System.Drawing.Color.AliceBlue;
            this.priceTest1.Controls.Add(this.priceResult);
            this.priceTest1.Controls.Add(this.listBox1);
            this.priceTest1.Controls.Add(this.pricingButton);
            this.priceTest1.Location = new System.Drawing.Point(8, 47);
            this.priceTest1.Name = "priceTest1";
            this.priceTest1.Padding = new System.Windows.Forms.Padding(3);
            this.priceTest1.Size = new System.Drawing.Size(1784, 743);
            this.priceTest1.TabIndex = 1;
            this.priceTest1.Text = "priceTest1";
            // 
            // priceResult
            // 
            this.priceResult.Location = new System.Drawing.Point(773, 186);
            this.priceResult.Name = "priceResult";
            this.priceResult.Size = new System.Drawing.Size(100, 46);
            this.priceResult.TabIndex = 2;
            // 
            // listBox1
            // 
            this.listBox1.FormattingEnabled = true;
            this.listBox1.ItemHeight = 32;
            this.listBox1.Location = new System.Drawing.Point(479, 168);
            this.listBox1.Name = "listBox1";
            this.listBox1.Size = new System.Drawing.Size(120, 68);
            this.listBox1.TabIndex = 1;
            // 
            // pricingButton
            // 
            this.pricingButton.Location = new System.Drawing.Point(66, 168);
            this.pricingButton.Name = "pricingButton";
            this.pricingButton.Size = new System.Drawing.Size(75, 42);
            this.pricingButton.TabIndex = 0;
            this.pricingButton.Text = "定";
            this.pricingButton.UseVisualStyleBackColor = true;
            // 
            // Form1
            // 
            this.AutoScaleDimensions = new System.Drawing.SizeF(13F, 24F);
            this.AutoScaleMode = System.Windows.Forms.AutoScaleMode.Font;
            this.ClientSize = new System.Drawing.Size(1774, 789);
            this.Controls.Add(this.Simulation);
            this.Name = "Form1";
            this.Text = "Form1";
            this.Simulation.ResumeLayout(false);
            this.priceTest2.ResumeLayout(false);
            this.pricingResult.ResumeLayout(false);
            this.pricingResult.PerformLayout();
            this.simulationPara.ResumeLayout(false);
            this.simulationPara.PerformLayout();
            this.msgGroupBox.ResumeLayout(false);
            this.hestonPara.ResumeLayout(false);
            this.hestonPara.PerformLayout();
            this.optionPara.ResumeLayout(false);
            this.optionPara.PerformLayout();
            this.tabPage1.ResumeLayout(false);
            this.tabPage1.PerformLayout();
            ((System.ComponentModel.ISupportInitialize)(this.chart1)).EndInit();
            this.priceTest1.ResumeLayout(false);
            this.priceTest1.PerformLayout();
            this.ResumeLayout(false);

        }

        #endregion
        private System.Windows.Forms.TabControl Simulation;
        private System.Windows.Forms.TabPage priceTest1;
        private System.Windows.Forms.TabPage tabPage1;
        private System.Windows.Forms.DataVisualization.Charting.Chart chart1;
        private System.Windows.Forms.RichTextBox richTextBox1;
        private System.Windows.Forms.Label label1;
        private System.Windows.Forms.Button button1;
        private System.Windows.Forms.TextBox textBox_rho_old;
        private System.Windows.Forms.Button pricingButton;
        private System.Windows.Forms.ListBox listBox1;
        private System.Windows.Forms.TextBox priceResult;
        private System.Windows.Forms.TabPage priceTest2;
        private System.Windows.Forms.GroupBox optionPara;
        private System.Windows.Forms.TextBox textBox_s0;
        private System.Windows.Forms.Label label3;
        private System.Windows.Forms.TextBox textBox_k;
        private System.Windows.Forms.Label label2;
        private System.Windows.Forms.TextBox textBox_var0;
        private System.Windows.Forms.Label label4;
        private System.Windows.Forms.Label label5;
        private System.Windows.Forms.TextBox textBox_rf;
        private System.Windows.Forms.TextBox textBox_T;
        private System.Windows.Forms.Label label6;
        private System.Windows.Forms.Button button2;
        private System.Windows.Forms.GroupBox hestonPara;
        private System.Windows.Forms.TextBox textBox_sigma;
        private System.Windows.Forms.Label label7;
        private System.Windows.Forms.Label label8;
        private System.Windows.Forms.TextBox textBox_theta;
        private System.Windows.Forms.Label label9;
        private System.Windows.Forms.Label label10;
        private System.Windows.Forms.TextBox textBox_kappa;
        private System.Windows.Forms.Label label11;
        private System.Windows.Forms.TextBox textBox_rho;
        private System.Windows.Forms.GroupBox msgGroupBox;
        private System.Windows.Forms.RichTextBox msgBox;
        private System.Windows.Forms.Button clear;
        private System.Windows.Forms.GroupBox simulationPara;
        private System.Windows.Forms.TextBox textBox_callPrice;
        private System.Windows.Forms.Label label13;
        private System.Windows.Forms.Label label12;
        private System.Windows.Forms.TextBox textBox_putPrice;
        private System.Windows.Forms.GroupBox pricingResult;
        private System.Windows.Forms.TextBox textBox_seed;
        private System.Windows.Forms.Label label15;
        private System.Windows.Forms.Label label14;
        private System.Windows.Forms.TextBox textBox_pathCnt;
    }
}

