function colourField(colour, screenRect)
    Screen('Preference', 'VisualDebugLevel', 1);
    myScreen = 0;
    window = Screen('OpenWindow',myScreen,colour,screenRect);
    Screen('Flip',window);
    
    % hold off starting the stimulation so I can have the shutter closed at
    % the beginning then open it just before I start stimulating
    KbWait;
    
    Screen('CloseAll');

    clear Screen;
end