
function [DH_entry] = gatherDHentry()

%A Basic file to take input from the user via a series of prompts

%close all

%Boolean for handling DH entry, 0  = false, 1 = true
%DH_entry = 0;


%Create a figure for the input ui
fig = uifigure('Name','DH Input','Position',[600 360 360 360]);
%Create a group for the selection buttons
bg = uibuttongroup(fig,'Position', [150 150 180 100]);
rb1 = uiradiobutton(bg,'Text','Enter DH Parameters','Position',[10 10 150 20], 'Value',0);
rb2 = uiradiobutton(bg,'Text','Enter Robot Parameters','Position',[10 30 180 20],'Value',0);

%Create a panel for the confirmation button
%pnl_con = uipanel(fig);
%Create a confirmaiton button on the panel
btn_con = uibutton(fig,'Text','Confirm','Position',[150 20 60 20],...
'ButtonPushedFcn', @(btn_con,event) ButtonPushed(fig, rb1,rb2));

%Start waiting to let the user adjust thier choice
uiwait(fig)
%Call the callback to find the value of DH_entry
ButtonPushed(fig,rb1,rb2)
%Close the input box
close(fig)

end

function DH_entry = ButtonPushed(fig, rb1, rb2) 

    %%Entering DH paramters directly
    if(rb1.Value == 1)
        DH_entry = 1;
        %resume after the confirm botton has been pressed
        uiresume(fig)
    %Enter Robot paramters
    elseif(rb2.Value == 1)
        DH_entry = 0;
        uiresume(fig)
    %Somehow if no buttons are selected
    else
        return
    end


end

