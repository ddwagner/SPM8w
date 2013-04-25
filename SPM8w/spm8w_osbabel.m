function PlatformString = spm8w_osbabel(UnixString)
% ==============================================================================
% SPM8w r5236
% Script driven batching for SPM8 with additional tools and support for 
% other commonly used analyses (roi, ppi, mixed).
% 
% Heatherton & Kelley Labs
% Last update: February 2013 - DDW
% Created: March, 2012 - DDW
% ==============================================================================
% FORMAT PlatformString = spm8w_osbabel(UnixString);
% 
% Simple function to translate unix system calls to windows specific calls.
% Feed it a string that would work on unix and it will return a string 
% for eval that will work on windows. 
% 
% This is useful as it obviates countless if/then statements in spm8w code.
% However, I've since realized that a lot of this can be done with matlab
% commands (delete, movefile, rmdir etc.) so I will slowly phase this out.
% ==============================================================================
% CHANGE LOG:
% =======1=========2=========3=========4=========5=========6=========7=========8


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Get OS (ispc, isunix, ismac)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ispc
    %Swap slashes
    PcString = strrep(UnixString, '/', filesep);
    %Check for cp, rm, mkdir, 
    [command, path] = strtok(PcString);
    switch(command)
        case 'cp'
            PcString = strrep(PcString, 'cp', 'copy');
        case 'mv'
            PcString = strrep(PcString, 'mv', 'move');
        case 'mkdir'
            %mkdir works on PC, nothing to do.
        case '!touch'
            PcString = ['!type nul>',path(2:end)];
        case 'rm'
            recurs = strtok(path); 
            if strcmp(recurs,'-r')
                PcString = strrep(PcString, 'rm -r', 'rmdir /S /Q');
            else  
                PcString = strrep(PcString, 'rm', 'del /Q');
            end
        otherwise
            if strcmp(command(2:3),':\') || strcmp(command(2:3),':/')
                %Nothing to do, we already fixed the slashes on line 24
            else
                error(['Unknown command for translation to PC: ',command]);
            end
    end
    PlatformString = PcString;
else
    PlatformString = UnixString;
end
