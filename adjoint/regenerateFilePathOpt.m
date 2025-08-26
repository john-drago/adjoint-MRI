function newFilePath = regenerateFilePathOpt( oldFilePath )

splitWord = "multiphoton";

tld = strsplit( fileparts(mfilename('fullpath')), splitWord );
tld = tld{ 1 };

filePathSplit = strsplit( oldFilePath, splitWord );
filePathSplitLower = strsplit( filePathSplit{ 2 }, [ "\", "/" ] );

newFilePath = fullfile( tld, splitWord, filePathSplitLower{ 2:end } );

end