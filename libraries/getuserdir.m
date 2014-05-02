function userDir = getuserdir
%GETUSERDIR   return the user home directory.
%   USERDIR = GETUSERDIR returns the user home directory using the registry
%   on windows systems and using Java on non windows systems as a string
%
%   Example:
%      getuserdir() revim gturns on windows
%           C:\Documents and Settings\MyName\Eigene Dateien
%
%   by Sven Probst, 10 Aug 2007 (Updated 10 Aug 2007) 
%   http://www.mathworks.de/matlabcentral/fileexchange/15885-get-user-home-directory

if ispc
    userDir = winqueryreg('HKEY_CURRENT_USER',...
        ['Software\Microsoft\Windows\CurrentVersion\' ...
         'Explorer\Shell Folders'],'Personal');
else
    userDir = char(java.lang.System.getProperty('user.home'));
end
