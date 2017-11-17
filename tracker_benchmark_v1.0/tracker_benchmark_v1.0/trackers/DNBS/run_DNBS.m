function results=run_DNBS(seq, res_path, bSaveImage)

close all;
disp('run-dnbs');
x=seq.init_rect(1)-1;%matlab to c
y=seq.init_rect(2)-1;
w=seq.init_rect(3);
h=seq.init_rect(4);

path = '.\results\';

%featureName kernelName param svmC svmBudgetSize searchRadius seed
%featureName: raw haar histogram
%kernelName: linear gaussian intersection chi2
%seed: default - 0
nz	= strcat('%0',num2str(seq.nz),'d');
style='.avi';
disp([seq.name style]);
tic
command = ['dnbsTracker.exe -box x=' num2str(x) ';y=' num2str(y) ';width=' num2str(w) ';height=' num2str(h) ' -from ' num2str(seq.startFrame) ' -to ' num2str(seq.endFrame)  ' -image ' seq.path  nz '.' seq.ext ' -set vers=9;nfgrd=3;ratio=0.5;iprod=0.75 -margin width=10;height=10 -savefile' ' ' seq.name ' -oavi ' seq.name style];
disp(command);
dos(command);
duration=toc;
results.type='rect';
results.fps=seq.len/duration;
disp([seq.name '_DNBS.txt']);
res=dlmread([seq.name '_DNBS.txt']);
 results.res = res;
 disp('txt³É¹¦');
 results.res(:,1:2) =results.res(:,1:2) + 1;%c to matlab
% 
% results.fps = dlmread([seq.name '_DNBS_FPS.txt']);

disp('results:');
disp(results)
disp(['fps: ' num2str(results.fps)]);
end

