function joinData()
fileNameBase = 'simpleMeasurement16Regions';
par_sets = {'1-5000',...
            '5001-8000',...
            '8001-10000',...
            '10001-15000',...
            '15001-20000',...
            '20001-30000',...
            };

for i = 1:length(par_sets)
    fileNames{i} =  [fileNameBase,'-',par_sets{i},'.mat'];
end
          
% fileNames = { 'simpleMeasurementAllTree-1-4000.mat', ...
%               'simpleMeasurementAllTree-4001-8000.mat', ...
%               'simpleMeasurementAllTree-8001-12000.mat', ...
%               'simpleMeasurementAllTree-12001-16000', ...
%               'simpleMeasurementAllTree-16001-20000', ...
%               'simpleMeasurementAllTree-20001-21000' };
          
% fileNames = { 'simpleMeasurementAllTreeNormalized-1-3000.mat', ...
%               'simpleMeasurementAllTreeNormalized-3001-6000.mat', ...
%               'simpleMeasurementAllTreeNormalized-6001-9000.mat', ...
%               'simpleMeasurementAllTreeNormalized-9001-12000', ...
%               'simpleMeasurementAllTreeNormalized-12001-15000', ...
%               'simpleMeasurementAllTreeNormalized-15001-18000', ...
%               'simpleMeasurementAllTreeNormalized-18001-19000', ...
%               'simpleMeasurementAllTreeNormalized-19001-20000', ...
%               'simpleMeasurementAllTreeNormalized-20001-30000'};

dataStruct.brainspanRandomResults = [];
dataStruct.brainspanResults = [];
dataStruct.human6PhysicalRandomResults = [];
dataStruct.human6PhysicalResults = [];
dataStruct.human6RandomResults  = [];
dataStruct.human6Results = [];

for i = 1:length(fileNames)
    currentDataStruct = load(fileNames{i} );
    dataStruct = joinAnotherData(dataStruct, currentDataStruct);
end
dataStruct.human_gene_info = currentDataStruct.human_gene_info;
dataStruct.developing_genes_info= currentDataStruct.developing_genes_info;

% assert(size(dataStruct.developing_genes_info,1) ==size(dataStruct.human6Results,1 ));
outputFile = [ fileNameBase, 'Par.mat'];
save( outputFile,'-struct','dataStruct');
end

function dataStruct = joinAnotherData(dataStruct, newDataStruct)
startIndex = newDataStruct.indices(1);
lastIndex = min(   newDataStruct.indices(2), size(newDataStruct.brainspanRandomResults,1) );
dataStruct.brainspanRandomResults = cat(1, dataStruct.brainspanRandomResults,newDataStruct.brainspanRandomResults(startIndex:lastIndex,:));
dataStruct.brainspanResults = cat(1, dataStruct.brainspanResults,newDataStruct.brainspanResults(startIndex:lastIndex,:));
lastIndex = min(   newDataStruct.indices(2), size(newDataStruct.human6PhysicalRandomResults,1) );
dataStruct.human6PhysicalRandomResults = cat(1, dataStruct.human6PhysicalRandomResults,newDataStruct.human6PhysicalRandomResults(startIndex:lastIndex ,:));
dataStruct.human6PhysicalResults = cat(1, dataStruct.human6PhysicalResults,newDataStruct.human6PhysicalResults(startIndex:lastIndex ,:));
lastIndex = min(   newDataStruct.indices(2), size(newDataStruct.human6RandomResults,1) );
dataStruct.human6RandomResults  = cat(1, dataStruct.human6RandomResults,newDataStruct.human6RandomResults(startIndex:lastIndex ,:));
dataStruct.human6Results = cat(1, dataStruct.human6Results,newDataStruct.human6Results(startIndex:lastIndex ,:));


end