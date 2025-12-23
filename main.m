%% Parallelization of Path Following Based Phase Unwrapping Algorithms

% 
% 1 - Input Images (Interferograms)
% 2 - Map Generation
% 3 - Phase Unwrapping Algorithms
%
% Dr. Emrah Onat
% 23.12.2025
% 

%% Main Code of InSAR

clear all
close all
clc

%% Display
disp('--- Algorithms ---');

fid = fopen( 'results.txt', 'wt' );
fprintf( fid, '%3s %13s %8s %8s %10s %10s %10s %10s\r\n','#i', '#Input','#Map', '#PUAlg','Duration', '#Residue', '#BranchCut', 'RMSE');

iteration = 0;
for i = 4:4
%     fprintf( fid, '%61s\r\n','------------------------------------------------------------');
    %% Input Images
    % 1 - P00 - ifsar.512x512
    % 9 - P00 - peaks.101x101
    % 16 -P00 - etna.1024x1024
    % 17 -P00 - etna.2000x5198
    
    numberofinputimage = i;
    [Inpimage, phaseimage, maskimage, corrimage, surfimage] = inputexamples(numberofinputimage, fid);
    
%     figure;
%     subplot(311);imagesc(phaseimage);title(['Input Phase Image, #Map = ' num2str(numberofinputimage)]);
%     subplot(312);imagesc(corrimage);title('Input Correlation Map');
%     subplot(313);mesh(surfimage);title('Groundtruth Unwrapped Map');

    for j = 5:5
        %% Guided Map Generation
        % 0 - No Map
        % 1 - Correlation Map
        % 2 - PseudoCorrelation Map
        % 3 - Phase Derivative Variance Map
        % 4 - Maximum Phase Variance Map
        % 5 - Amplitude Variance Map
        
        numberofQualAlgo = j;
        [MAPtype, qualmap, average_val] = QualMapGen(numberofQualAlgo, phaseimage, corrimage, maskimage, fid);
%         figure;imagesc(qualmap);title(['Guide Map, Average Value = ' num2str(average_val) ', Map Number = ' num2str(numberofQualAlgo)]);

        for k = 1:1
            %% Phase Unwrapping Algorithms

            % 2  - Goldstein C/C++
            % 4  - Quality Guided C/C++
            % 6  - Flynn C/C++
			% 19 - HBP    (My Algo: "Phase unwrapping via hierarchical and balanced residue partitioning")
            
            if k==7 && (i==5 || i== 6 || i==7 || i==8 || i==1 || i==11)
                continue;
            end
            if k==6 && (i==10)
                continue;
            end
            if k==3 
                continue;
            end

            numberofPUAlgo = k;
            tic;
            [PUAlg, resmap, BCmap, unwrappedmap] = PUalgorithms(numberofPUAlgo, phaseimage, maskimage, qualmap);
            duration = toc;
            resnumber = length(find(resmap));
            BClength = sum(sum(BCmap));
            unwrappedmap = unwrappedmap-min(min(unwrappedmap)); unwrappedmap = unwrappedmap/max(max(unwrappedmap)); 
            surfimage = surfimage-min(min(surfimage)); surfimage = surfimage/max(max(surfimage)); 
            rmse_uW = rms(rms(surfimage-unwrappedmap));
            
            figure
            subplot(321);imagesc(phaseimage);title(['Input Phase Image, #Input = ' Inpimage]);
            subplot(322);imagesc(qualmap);title(['Guide Map, Avg Val = ' num2str(average_val) ', #Map = ' MAPtype]);
            subplot(323);imagesc(resmap);title(['Residue Map, #Res = ' num2str(resnumber)]);
            subplot(324);imagesc(BCmap);title(['Branch-Cut Map, length = ' num2str(BClength)]);
            subplot(325);imagesc(unwrappedmap);title(['Unwrapped Map, RMSE = ' num2str(rmse_uW)]);
            subplot(326);mesh(unwrappedmap);title(['Unwrapped Map, PU Algo = ' PUAlg]);
            
            iteration = iteration +1;
            A = [iteration; i; j; k; resnumber; BClength; rmse_uW; duration];

            fprintf( fid, '%3d %13s %8s %8s %10f %10d %10d %10f\r\n', iteration, Inpimage, MAPtype, PUAlg, duration, resnumber, BClength, rmse_uW);
%             fprintf( fid, '%5s \r\n', PUAlg);
        end
    end
end

fclose(fid);
open('results.txt');