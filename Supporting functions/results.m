function dummyOutput = results(tr)
    % Give user options to save or view results (or exit)
    msg = 'View result';
    while strcmp(msg,'View result') || strcmp(msg,'Export STL')
        msg = questdlg(["Topology optimization successful!", "Select an option from the menu below."],...
            "Results", "View result", "Export STL", "Finish", "View result");

        % Select options based on user input (proceed only if "Proceed" selected)
        switch msg
            case {'Finish', ''}

                % As a precaution, check whether the user intended to exit
                exitCfm = "Help";
                while strcmp(exitCfm,"Help")
                    exitCfm = questdlg(["Are you sure you want to finish?", "Your current settings " + ...
                        "and results will be lost."], "Ready to leave?", "Yes, I'm done", "No, go back", ...
                        "Help", "Yes, I'm done");

                    % Select next step based on user input
                    switch exitCfm

                        case "Help"

                            % Display help text
                            helpBox = msgbox("To keep your workspace clean, the optio program will " + ...
                                "automatically clear your workspace variables on exit. We do not " + ...
                                "currently support any exporting of workspace variables. If you are " + ...
                                "unsure of whether you need these results, we recommend that you " + ...
                                "save a copy of the STL as a precaution.", "Help", 'help');
                            uiwait(helpBox);

                        case {"No, go back", ''}

                            % Reassign msg to redirect the user back to
                            % the results page
                            msg = 'View result';
                    end
                end

            case 'View result'

                % Create new figure window and plot mesh
                result = figure('Name','Optimization result','NumberTitle','off');
                trimesh(tr);
                axis off;
                axis image;

                % Select 3D rotation tool by default and wait until user closes
                % the window to return to the results window
                rotate3d on;
                uiwait(result);

            case 'Export STL'

                % Open a dialog box for filesaving
                [file,path] = uiputfile(date + "_optimization-result.stl");

                % Proceed only if user selected a path and filename
                if ~isequal(file,0) && ~isequal(path,0)

                    % Export the mesh as an STL to the appropriate path, with
                    % the correct filename
                    stlwrite(tr,fullfile(path,file));

                    % Notify sucessful save
                    savemsg('File saved!');
                    continue
                end

                % Notify aborted save
                savemsg('Aborted save.');
        end
    end
    
    dummyOutput = NaN;
end