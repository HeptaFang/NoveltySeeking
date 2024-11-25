classdef DataGUI_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                 matlab.ui.Figure
        GridLayout               matlab.ui.container.GridLayout
        KernelButtonGroup        matlab.ui.container.ButtonGroup
        Kernel1Button            matlab.ui.control.RadioButton
        Kernel2Button            matlab.ui.control.RadioButton
        Kernel3Button            matlab.ui.control.RadioButton
        TestOutput               matlab.ui.control.Label
        TrialSlider              matlab.ui.control.Slider
        TrialSliderLabel         matlab.ui.control.Label
        Neuron1SliderLabel       matlab.ui.control.Label
        Neuron1Slider            matlab.ui.control.Slider
        Neuron2SliderLabel       matlab.ui.control.Label
        Neuron2Slider            matlab.ui.control.Slider
        SessionTypeDropDown      matlab.ui.control.DropDown
        SessionDropDownLabel     matlab.ui.control.Label
        TrialSpinner             matlab.ui.control.Spinner
        Neuron1Spinner           matlab.ui.control.Spinner
        Neuron2Spinner           matlab.ui.control.Spinner
        ConnectionDispLabel      matlab.ui.control.Label
        PrePostButtonGroup       matlab.ui.container.ButtonGroup
        PreButton                matlab.ui.control.RadioButton
        PostButton               matlab.ui.control.RadioButton
        SessionStageButtonGroup  matlab.ui.container.ButtonGroup
        RestingStagesLabel       matlab.ui.control.Label
        TaskStagesLabel          matlab.ui.control.Label
        RestCloseButton          matlab.ui.control.RadioButton
        RestOpenButton           matlab.ui.control.RadioButton
        InfoRespButton           matlab.ui.control.RadioButton
        InfoButton               matlab.ui.control.RadioButton
        InfoAntiButton           matlab.ui.control.RadioButton
        DecisionButton           matlab.ui.control.RadioButton
        LoadButton               matlab.ui.control.Button
        LoadingStatusLabel       matlab.ui.control.Label
        SessionIdxDropDown       matlab.ui.control.DropDown
        ConnectionAxes           matlab.ui.control.UIAxes
        RawAxes                  matlab.ui.control.UIAxes
    end

    properties (Access = private)
        % constants
        all_session_types = {'Muscimol', 'Saline', 'SimRec'};

        % file path constants/variables
        root_path = '\\storage1.ris.wustl.edu\ilyamonosov\Active\Tianhong\'
        kernel = 'Delta';
        reg = 'L2=3e3';
        epoch = '2500';

        % environment options
        cmap = NaN;
        colorlim = 2;
    
        % status variables
        file_found = false;
        file_loaded = false;
        conn_plotted = false;

        % option variables
        session_type = '';
        session_idx = '-';
        session_prepost = '';
        session_stage = '';

        current_kernel = 1;
        current_trial = 1; % update this.

        cursor_x = 0;
        cursor_y = 0;

        % data
        J_mat = NaN;
        N = 0;
        borders = NaN;
        cell_id = NaN;
        cell_area = NaN;
        sort_idx = NaN;
        rasters = NaN;
        spikes = NaN;
        n_trial = 0;
    end
    
    methods (Access = private)
        function SyncVariables(app)
            app.session_type = app.SessionTypeDropDown.Value;
            app.session_idx = app.SessionIdxDropDown.Value;
            app.session_prepost = app.PrePostButtonGroup.SelectedObject.Text;
            app.session_stage = app.SessionStageButtonGroup.SelectedObject.Text;
        end

        function UpdateSessionIdxSelection(app)
            session_idx_list = [];
            switch app.session_type
                case 'Muscimol'
                    session_idx_list = 1:10;
                case 'Saline'
                    session_idx_list = 1:5;
                case 'SimRec'
                    session_idx_list = 1:5;
            end

            session_idx_cell = {'-'};
            for i = session_idx_list
                session_idx_cell{end+1} = int2str(i);
            end

            app.SessionIdxDropDown.Items = session_idx_cell;
            app.SessionIdxDropDown.Value = '-';
            app.session_idx = '-';
            app.UpdateSessionStageSelection();
        end

        function UpdateSessionStageSelection(app)
            if strcmp(app.session_idx, '-')
                app.PrePostButtonGroup.Enable = "off";
                app.SessionStageButtonGroup.Enable = "off";
            else
                app.PrePostButtonGroup.Enable = "on";
                app.SessionStageButtonGroup.Enable = "on";
            end
        end

        function UpdateLoadButton(app)
            if strcmp(app.session_idx, '-')
                app.LoadButton.Enable = "off";
                app.LoadingStatusLabel.Text = "Select a Session.";
            else
                % find file
                app.file_found = false;
                if strcmp(app.session_prepost, 'Pre')
                    session_stage_full = [app.session_type, 'Pre', app.session_stage, '_full'];
                else
                    session_stage_full = [app.session_type, 'Post', app.session_stage, '_cortex'];
                end

                file_path = [app.root_path, 'GLM_model/', session_stage_full,...
                    '/GLM_', session_stage_full, '_', app.session_idx, '_',...
                    app.kernel, '_0_', app.reg, '_', app.epoch, '.mat'];
                app.TestOutput.Text = file_path;
                app.file_found = isfile(file_path);

                if app.file_found
                    app.LoadButton.Enable = "on";
                    app.LoadingStatusLabel.Text = "Ready to Load.";
                else
                    app.LoadButton.Enable = "off";
                    app.LoadingStatusLabel.Text = "File Not Found.";
                end
            end
        end

        function PlotConnection(app)
            image = app.J_mat(:, :, app.current_kernel);
            imagesc(app.ConnectionAxes, image, "ButtonDownFcn", createCallbackFcn(app, @ConnectionAxesButtonDown, true));
            clim(app.ConnectionAxes, [-app.colorlim, app.colorlim]);
            colormap(app.ConnectionAxes, app.cmap);
            axis(app.ConnectionAxes, "square");
            xlim(app.ConnectionAxes, [0.5, app.N+0.5]);
            ylim(app.ConnectionAxes, [0.5, app.N+0.5]);
            app.UpdateCursor();

            % border lines
            for b=app.borders
                line(app.ConnectionAxes, [0.5, app.N+0.5], [b, b], 'Color', 'k');
                line(app.ConnectionAxes, [b, b], [0.5, app.N+0.5], 'Color', 'k');
            end

            app.conn_plotted = true;
            app.Neuron1Slider.Enable = "on";
            app.Neuron2Slider.Enable = "on";
            app.Neuron1Spinner.Enable = "on";
            app.Neuron2Spinner.Enable = "on";
            if app.n_trial > 1
                app.TrialSlider.Enable = "on";
                app.TrialSpinner.Enable = "on";
            else
                app.TrialSlider.Enable = "off";
                app.TrialSpinner.Enable = "off";
            end
            app.KernelButtonGroup.Enable = "on";
            app.ConnectionAxes.ButtonDownFcn = createCallbackFcn(app, @ConnectionAxesButtonDown, true);
        end

        function PlotRaster(app)
            cla(app.RawAxes);
            trial_activity = app.rasters{app.current_trial};
            hold(app.RawAxes, "on");
            plot(app.RawAxes, 2 + trial_activity(app.cursor_x, :));
            plot(app.RawAxes, 0 + trial_activity(app.cursor_y, :));
            hold(app.RawAxes, "off");
            ylim(app.RawAxes, [-0.5, 3.5]);
            % z = zoom(app.RawAxes);
            % setAxesZoomConstraint(z, app.RawAxes, 'x');
            % z.Motion = 'horizontal';
            zoom(app.RawAxes, 'xon');
            zoom(app.RawAxes, 'off');
        end

        function UpdateCursor(app)
            app.cursor_x = app.Neuron1Slider.Value;
            app.cursor_y = app.Neuron2Slider.Value;    
            hold(app.ConnectionAxes, "on");
            % xmin = app.cursor_x - 0.5;
            % ymin = app.cursor_y - 0.5;
            % rectangle(app.ConnectionAxes, 'position',[xmin ymin 1 1],'LineWidth',0.1)
            rectangle(app.ConnectionAxes, 'position',[app.cursor_x-0.5 0.5 1 app.N+0.5],'LineWidth',0.1, 'EdgeColor', 'r')
            rectangle(app.ConnectionAxes, 'position',[0.5 app.cursor_y-0.5 app.N+0.5 1],'LineWidth',0.1, 'EdgeColor', 'r')
            hold(app.ConnectionAxes, "off");
            app.UpdateDisplay();
            app.PlotRaster();
        end

        function UpdateDisplay(app)
            connection_string = sprintf("%s %s to %s %s, kernel %d, J = %0.3f\naveJ = %d", ...
                app.cell_area{app.cursor_x}, app.cell_id{app.cursor_x}, ...
                app.cell_area{app.cursor_y}, app.cell_id{app.cursor_y}, ...
                app.current_kernel, app.J_mat(app.cursor_y, app.cursor_x, app.current_kernel), ...
                mean(app.J_mat, "all"));
            app.ConnectionDispLabel.Text = connection_string;
        end

        function SyncSliderToSpinner(app)
            app.Neuron1Spinner.Value = app.Neuron1Slider.Value;
            app.Neuron2Spinner.Value = app.Neuron2Slider.Value;
            app.TrialSpinner.Value = app.TrialSlider.Value;
        end

        function SyncSpinnerToSlider(app)
            app.Neuron1Slider.Value = app.Neuron1Spinner.Value;
            app.Neuron2Slider.Value = app.Neuron2Spinner.Value;
            app.TrialSlider.Value = app.TrialSpinner.Value;
        end
    end
    

    % Callbacks that handle component events
    methods (Access = private)

        % Code that executes after component creation
        function startupFcn(app)
            % Initialize session types
            app.SessionTypeDropDown.Items = app.all_session_types;
            app.SessionTypeDropDown.Value = app.all_session_types{1};
            app.session_type = app.all_session_types{1};

            % load colormap
            app.cmap=brewermap(256,'*RdBu');

            % init axes
            axis(app.ConnectionAxes, "square");

            app.UpdateSessionIdxSelection();
            app.UpdateLoadButton();
        end

        % Value changed function: SessionTypeDropDown
        function SessionTypeDropDownValueChanged(app, event)
            app.SyncVariables();
            app.UpdateSessionIdxSelection();
            app.UpdateLoadButton();
        end

        % Value changed function: SessionIdxDropDown
        function SessionIdxDropDownValueChanged(app, event)
            app.SyncVariables();
            app.UpdateSessionStageSelection();
            app.UpdateLoadButton();
        end

        % Selection changed function: PrePostButtonGroup, 
        % ...and 1 other component
        function SessionStageButtonGroupSelectionChanged(app, event)
            app.SyncVariables();
            app.UpdateLoadButton();
        end

        % Button pushed function: LoadButton
        function LoadButtonPushed(app, event)
            if strcmp(app.session_prepost, 'Pre')
                session_stage_full = [app.session_type, 'Pre', app.session_stage, '_full'];
            else
                session_stage_full = [app.session_type, 'Post', app.session_stage, '_cortex'];
            end

            % load model
            file_path = [app.root_path, 'GLM_model/', session_stage_full,...
                    '/GLM_', session_stage_full, '_', app.session_idx, '_',...
                    app.kernel, '_0_', app.reg, '_', app.epoch, '.mat'];
            load(file_path, "model_par", "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel", "kernel_len", "N")
            app.N = N;
            app.J_mat = zeros(N, N, n_conn_kernel);
            
            % load metadata
            load([app.root_path, 'GLM_data/', session_stage_full,'/borders_', session_stage_full, '_', ...
                app.session_idx,'.mat'], "borders");
            app.borders = borders;
            load([app.root_path, 'GLM_data/', session_stage_full,'/sortidx_', session_stage_full, '_', ...
                app.session_idx, '_', app.kernel, '.mat'], "sort_idx");
            app.sort_idx = sort_idx;
            load([app.root_path, 'GLM_data/', session_stage_full,'/raster_',session_stage_full, '_', ...
                app.session_idx,'_0.mat'], "cell_area", "cell_id", "rasters", "spikes", "n_trial");
            app.cell_area = cell_area;
            app.cell_id = cell_id;
            app.rasters = rasters;
            app.spikes = spikes;
            app.n_trial = n_trial;

            file_path = [app.root_path, 'GLM_model/', session_stage_full,...
                    '/GLM_', session_stage_full, '_', app.session_idx, '_',...
                    app.kernel, '_0_', app.reg, '_', app.epoch, '.mat'];
            load(file_path, "model_par", "PS_kernels", "conn_kernels", "n_PS_kernel", "n_conn_kernel", "kernel_len", "N")

            for i = 1:n_conn_kernel
                app.J_mat(:, :, i) = model_par(:, (app.N*(i-1) + n_PS_kernel + 2):(app.N*i + n_PS_kernel + 1));
                app.J_mat(:, :, i) = app.J_mat(app.sort_idx, app.sort_idx, i);
            end

            % sort
            app.cell_area = app.cell_area(app.sort_idx);
            app.cell_id = app.cell_id(app.sort_idx);
            app.spikes = app.spikes(app.sort_idx, :);
            for i = 1:app.n_trial
                app.rasters{i} = app.rasters{i}(app.sort_idx, :);
            end

            app.Neuron1Slider.Value = 1;
            app.Neuron2Slider.Value = 1;
            app.Neuron1Slider.Limits = [1, app.N];
            app.Neuron2Slider.Limits = [1, app.N];
            app.Neuron1Slider.MinorTicks = 1:app.N;
            app.Neuron2Slider.MinorTicks = 1:app.N;
            app.TrialSlider.Value = 1;
            if app.n_trial > 1
                app.TrialSlider.Enable = "on";
                app.TrialSpinner.Enable = "on";
                app.TrialSlider.Limits = [1, app.n_trial];
                app.TrialSlider.MinorTicks = 1:app.n_trial;
            else
                app.TrialSlider.Enable = "off";
                app.TrialSpinner.Enable = "off";
                app.TrialSlider.Limits = [0, 1];
                app.TrialSlider.MinorTicks = 0:1;
            end

            % app.current_kernel = 1;
            app.current_trial = 1;
            app.cursor_x = 1;
            app.cursor_y = 1;

            app.SyncSliderToSpinner();
            app.PlotConnection();
        end

        % Callback function: Neuron1Slider, Neuron1Slider, Neuron2Slider, 
        % ...and 3 other components
        function SliderValueChanged(app, event)
            app.Neuron1Slider.Value = round(app.Neuron1Slider.Value);
            app.Neuron2Slider.Value = round(app.Neuron2Slider.Value);
            app.TrialSlider.Value = round(app.TrialSlider.Value);
            app.cursor_x = app.Neuron1Slider.Value;
            app.cursor_y = app.Neuron2Slider.Value;
            app.current_trial = app.TrialSlider.Value;
            app.SyncSliderToSpinner();
            app.PlotConnection();
        end

        % Value changed function: Neuron1Spinner, Neuron2Spinner, 
        % ...and 1 other component
        function SpinnerValueChanged(app, event)
            app.Neuron1Spinner.Value = round(app.Neuron1Spinner.Value);
            app.Neuron2Spinner.Value = round(app.Neuron2Spinner.Value);
            app.TrialSpinner.Value = round(app.TrialSpinner.Value);
            app.cursor_x = app.Neuron1Spinner.Value;
            app.cursor_y = app.Neuron2Spinner.Value;
            app.current_trial = app.TrialSpinner.Value;
            app.SyncSpinnerToSlider();
            app.PlotConnection();
        end

        % Callback function
        function ConnectionAxesButtonDown(app, event)
            % fprintf("click, %0.3f, %0.3f\n", event.IntersectionPoint(1), event.IntersectionPoint(2));
            converted_x = round(event.IntersectionPoint(1));
            converted_y = round(event.IntersectionPoint(2));
            app.Neuron1Spinner.Value = converted_x;
            app.Neuron2Spinner.Value = converted_y;
            app.cursor_x = app.Neuron1Spinner.Value;
            app.cursor_y = app.Neuron2Spinner.Value;
            app.SyncSpinnerToSlider();
            app.PlotConnection();
        end

        % Selection changed function: KernelButtonGroup
        function KernelButtonGroupSelectionChanged(app, event)
            selectedButton = app.KernelButtonGroup.SelectedObject;
            switch selectedButton.Text
                case 'Kernel 1'
                    app.current_kernel = 1;
                case 'Kernel 2'
                    app.current_kernel = 2;
                case 'Kernel 3'
                    app.current_kernel = 3;
            end

            app.PlotConnection();
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.Position = [100 100 942 706];
            app.UIFigure.Name = 'MATLAB App';

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {'6x', '6x', '6x', '6x', '6x', '2x', '3x', '30x', '6x'};
            app.GridLayout.RowHeight = {'4.14x', '6x', '1x', 37, '1x', 40, 49, 22, 18, 40, 23};
            app.GridLayout.ColumnSpacing = 8.5;
            app.GridLayout.RowSpacing = 3;
            app.GridLayout.Padding = [8.5 3 8.5 3];

            % Create RawAxes
            app.RawAxes = uiaxes(app.GridLayout);
            title(app.RawAxes, 'Raster')
            xlabel(app.RawAxes, 'Time (ms)')
            zlabel(app.RawAxes, 'Z')
            app.RawAxes.YLim = [-0.5 3.5];
            app.RawAxes.YTick = [0.5 2.5];
            app.RawAxes.YTickLabel = {'Neuron 2'; 'Neuron 1'};
            app.RawAxes.Layout.Row = [1 2];
            app.RawAxes.Layout.Column = [6 9];

            % Create ConnectionAxes
            app.ConnectionAxes = uiaxes(app.GridLayout);
            title(app.ConnectionAxes, 'Connection Matrix')
            app.ConnectionAxes.PlotBoxAspectRatio = [1 1 1];
            app.ConnectionAxes.XTick = [];
            app.ConnectionAxes.YTick = [];
            app.ConnectionAxes.Layout.Row = [1 5];
            app.ConnectionAxes.Layout.Column = [1 5];

            % Create SessionIdxDropDown
            app.SessionIdxDropDown = uidropdown(app.GridLayout);
            app.SessionIdxDropDown.Items = {'1'};
            app.SessionIdxDropDown.ValueChangedFcn = createCallbackFcn(app, @SessionIdxDropDownValueChanged, true);
            app.SessionIdxDropDown.Layout.Row = 6;
            app.SessionIdxDropDown.Layout.Column = 5;
            app.SessionIdxDropDown.Value = '1';

            % Create LoadingStatusLabel
            app.LoadingStatusLabel = uilabel(app.GridLayout);
            app.LoadingStatusLabel.Layout.Row = 11;
            app.LoadingStatusLabel.Layout.Column = [2 3];
            app.LoadingStatusLabel.Text = 'Select a Session.';

            % Create LoadButton
            app.LoadButton = uibutton(app.GridLayout, 'push');
            app.LoadButton.ButtonPushedFcn = createCallbackFcn(app, @LoadButtonPushed, true);
            app.LoadButton.Enable = 'off';
            app.LoadButton.Layout.Row = 11;
            app.LoadButton.Layout.Column = 4;
            app.LoadButton.Text = 'Load';

            % Create SessionStageButtonGroup
            app.SessionStageButtonGroup = uibuttongroup(app.GridLayout);
            app.SessionStageButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @SessionStageButtonGroupSelectionChanged, true);
            app.SessionStageButtonGroup.Enable = 'off';
            app.SessionStageButtonGroup.Title = 'Stage';
            app.SessionStageButtonGroup.Layout.Row = [8 10];
            app.SessionStageButtonGroup.Layout.Column = [1 5];

            % Create DecisionButton
            app.DecisionButton = uiradiobutton(app.SessionStageButtonGroup);
            app.DecisionButton.Text = 'Decision';
            app.DecisionButton.Position = [11 19 68 22];
            app.DecisionButton.Value = true;

            % Create InfoAntiButton
            app.InfoAntiButton = uiradiobutton(app.SessionStageButtonGroup);
            app.InfoAntiButton.Text = 'InfoAnti';
            app.InfoAntiButton.Position = [11 -2 65 22];

            % Create InfoButton
            app.InfoButton = uiradiobutton(app.SessionStageButtonGroup);
            app.InfoButton.Text = 'Info';
            app.InfoButton.Position = [81 19 65 22];

            % Create InfoRespButton
            app.InfoRespButton = uiradiobutton(app.SessionStageButtonGroup);
            app.InfoRespButton.Text = 'InfoResp';
            app.InfoRespButton.Position = [81 -2 70 22];

            % Create RestOpenButton
            app.RestOpenButton = uiradiobutton(app.SessionStageButtonGroup);
            app.RestOpenButton.Text = 'RestOpen';
            app.RestOpenButton.Position = [182 -2 76 22];

            % Create RestCloseButton
            app.RestCloseButton = uiradiobutton(app.SessionStageButtonGroup);
            app.RestCloseButton.Text = 'RestClose';
            app.RestCloseButton.Position = [182 19 77 22];

            % Create TaskStagesLabel
            app.TaskStagesLabel = uilabel(app.SessionStageButtonGroup);
            app.TaskStagesLabel.Position = [13 40 70 22];
            app.TaskStagesLabel.Text = 'Task Stages';

            % Create RestingStagesLabel
            app.RestingStagesLabel = uilabel(app.SessionStageButtonGroup);
            app.RestingStagesLabel.Position = [182 40 86 22];
            app.RestingStagesLabel.Text = 'Resting Stages';

            % Create PrePostButtonGroup
            app.PrePostButtonGroup = uibuttongroup(app.GridLayout);
            app.PrePostButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @SessionStageButtonGroupSelectionChanged, true);
            app.PrePostButtonGroup.Enable = 'off';
            app.PrePostButtonGroup.Title = 'Pre/Post';
            app.PrePostButtonGroup.Layout.Row = 7;
            app.PrePostButtonGroup.Layout.Column = [1 5];

            % Create PostButton
            app.PostButton = uiradiobutton(app.PrePostButtonGroup);
            app.PostButton.Text = 'Post';
            app.PostButton.Position = [161 3 65 22];

            % Create PreButton
            app.PreButton = uiradiobutton(app.PrePostButtonGroup);
            app.PreButton.Text = 'Pre';
            app.PreButton.Position = [46 3 58 22];
            app.PreButton.Value = true;

            % Create ConnectionDispLabel
            app.ConnectionDispLabel = uilabel(app.GridLayout);
            app.ConnectionDispLabel.HorizontalAlignment = 'center';
            app.ConnectionDispLabel.FontSize = 16;
            app.ConnectionDispLabel.FontWeight = 'bold';
            app.ConnectionDispLabel.Layout.Row = 3;
            app.ConnectionDispLabel.Layout.Column = [6 9];
            app.ConnectionDispLabel.Text = '';

            % Create Neuron2Spinner
            app.Neuron2Spinner = uispinner(app.GridLayout);
            app.Neuron2Spinner.ValueChangedFcn = createCallbackFcn(app, @SpinnerValueChanged, true);
            app.Neuron2Spinner.Enable = 'off';
            app.Neuron2Spinner.Layout.Row = 10;
            app.Neuron2Spinner.Layout.Column = 9;

            % Create Neuron1Spinner
            app.Neuron1Spinner = uispinner(app.GridLayout);
            app.Neuron1Spinner.ValueChangedFcn = createCallbackFcn(app, @SpinnerValueChanged, true);
            app.Neuron1Spinner.Enable = 'off';
            app.Neuron1Spinner.Layout.Row = [8 9];
            app.Neuron1Spinner.Layout.Column = 9;

            % Create TrialSpinner
            app.TrialSpinner = uispinner(app.GridLayout);
            app.TrialSpinner.ValueChangedFcn = createCallbackFcn(app, @SpinnerValueChanged, true);
            app.TrialSpinner.Enable = 'off';
            app.TrialSpinner.Layout.Row = 7;
            app.TrialSpinner.Layout.Column = 9;

            % Create SessionDropDownLabel
            app.SessionDropDownLabel = uilabel(app.GridLayout);
            app.SessionDropDownLabel.HorizontalAlignment = 'center';
            app.SessionDropDownLabel.Layout.Row = 6;
            app.SessionDropDownLabel.Layout.Column = 1;
            app.SessionDropDownLabel.Text = 'Session';

            % Create SessionTypeDropDown
            app.SessionTypeDropDown = uidropdown(app.GridLayout);
            app.SessionTypeDropDown.Items = {'Muscimol', 'Saline', 'SimRec'};
            app.SessionTypeDropDown.ValueChangedFcn = createCallbackFcn(app, @SessionTypeDropDownValueChanged, true);
            app.SessionTypeDropDown.Layout.Row = 6;
            app.SessionTypeDropDown.Layout.Column = [2 4];
            app.SessionTypeDropDown.Value = 'Muscimol';

            % Create Neuron2Slider
            app.Neuron2Slider = uislider(app.GridLayout);
            app.Neuron2Slider.Limits = [1 10];
            app.Neuron2Slider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Neuron2Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Neuron2Slider.Enable = 'off';
            app.Neuron2Slider.Layout.Row = 10;
            app.Neuron2Slider.Layout.Column = 8;
            app.Neuron2Slider.Value = 1;

            % Create Neuron2SliderLabel
            app.Neuron2SliderLabel = uilabel(app.GridLayout);
            app.Neuron2SliderLabel.HorizontalAlignment = 'center';
            app.Neuron2SliderLabel.Layout.Row = 10;
            app.Neuron2SliderLabel.Layout.Column = [6 7];
            app.Neuron2SliderLabel.Text = 'Neuron 2';

            % Create Neuron1Slider
            app.Neuron1Slider = uislider(app.GridLayout);
            app.Neuron1Slider.Limits = [1 10];
            app.Neuron1Slider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Neuron1Slider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.Neuron1Slider.Enable = 'off';
            app.Neuron1Slider.Layout.Row = [8 9];
            app.Neuron1Slider.Layout.Column = 8;
            app.Neuron1Slider.Value = 1;

            % Create Neuron1SliderLabel
            app.Neuron1SliderLabel = uilabel(app.GridLayout);
            app.Neuron1SliderLabel.HorizontalAlignment = 'center';
            app.Neuron1SliderLabel.Layout.Row = [8 9];
            app.Neuron1SliderLabel.Layout.Column = [6 7];
            app.Neuron1SliderLabel.Text = 'Neuron 1';

            % Create TrialSliderLabel
            app.TrialSliderLabel = uilabel(app.GridLayout);
            app.TrialSliderLabel.HorizontalAlignment = 'center';
            app.TrialSliderLabel.Layout.Row = 7;
            app.TrialSliderLabel.Layout.Column = [6 7];
            app.TrialSliderLabel.Text = 'Trial';

            % Create TrialSlider
            app.TrialSlider = uislider(app.GridLayout);
            app.TrialSlider.Limits = [1 10];
            app.TrialSlider.ValueChangedFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.TrialSlider.ValueChangingFcn = createCallbackFcn(app, @SliderValueChanged, true);
            app.TrialSlider.Enable = 'off';
            app.TrialSlider.Layout.Row = 7;
            app.TrialSlider.Layout.Column = 8;
            app.TrialSlider.Value = 1;

            % Create TestOutput
            app.TestOutput = uilabel(app.GridLayout);
            app.TestOutput.WordWrap = 'on';
            app.TestOutput.Enable = 'off';
            app.TestOutput.Visible = 'off';
            app.TestOutput.Layout.Row = 5;
            app.TestOutput.Layout.Column = 9;
            app.TestOutput.Text = 'Path';

            % Create KernelButtonGroup
            app.KernelButtonGroup = uibuttongroup(app.GridLayout);
            app.KernelButtonGroup.SelectionChangedFcn = createCallbackFcn(app, @KernelButtonGroupSelectionChanged, true);
            app.KernelButtonGroup.Enable = 'off';
            app.KernelButtonGroup.Title = 'Kernel';
            app.KernelButtonGroup.Layout.Row = [5 6];
            app.KernelButtonGroup.Layout.Column = 8;

            % Create Kernel3Button
            app.Kernel3Button = uiradiobutton(app.KernelButtonGroup);
            app.Kernel3Button.Text = 'Kernel 3';
            app.Kernel3Button.Position = [259 23 67 22];

            % Create Kernel2Button
            app.Kernel2Button = uiradiobutton(app.KernelButtonGroup);
            app.Kernel2Button.Text = 'Kernel 2';
            app.Kernel2Button.Position = [148 23 67 22];

            % Create Kernel1Button
            app.Kernel1Button = uiradiobutton(app.KernelButtonGroup);
            app.Kernel1Button.Text = 'Kernel 1';
            app.Kernel1Button.Position = [37 23 67 22];
            app.Kernel1Button.Value = true;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = DataGUI_exported

            % Create UIFigure and components
            createComponents(app)

            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            % Execute the startup function
            runStartupFcn(app, @startupFcn)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end