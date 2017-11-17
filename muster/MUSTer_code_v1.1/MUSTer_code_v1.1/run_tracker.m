 
tic;
addpath('mexopencv');

addpath('ICF');
base_path = './';
res_path = 'Results/';

name = 'Shaking';

% video_path = [base_path name '/'];
video_path = 'C:\Users\ww\Desktop\visualtracking\sequences\Matrix\';

[ source.img_files, pos, target_sz, ground_truth, source.video_path]...
 = load_video_info(video_path);
source.n_frames = numel(source.img_files);
rect_init = [pos, target_sz];
% 
bboxes = MUSTer_tracking(source, rect_init);

dlmwrite([res_path name '.txt'], bboxes);
t=toc;
disp(source.n_frames);
fps=source.n_frames/t;
disp(fps);

writerObj=VideoWriter('jogging_muster.avi');
open(writerObj);
boxes=dlmread([res_path name '.txt']);

for i=1:source.n_frames
    img=imread([source.video_path source.img_files{i}]);
    imshow(img);
    hold on
    initstate=[boxes(i,1) boxes(i,2) boxes(i,3) boxes(i,4)];
     rectangle('Position',initstate,'LineWidth',4,'EdgeColor','r');
     frame=getframe(gcf);
     writeVideo(writerObj,frame.cdata);
end
close(writerObj);

