function spm_bet(Pt1)

%% BET brain and remove lose voxels from parenchyma
% spm_bet(Pt1)
% Based on spm_clean_brain
% Load native Pt1 space image(s) to create brain
% mask. The brain must have been segmented first in SPM, or segments
% created with SPM naming convention c1image, c2image ...
% Either pass a list of T1 images to be betted or you will be prompteed to
% select files.
% Amended to include 5 or 6 segments
% Mask t1 map with brain mask
% Output files according to the the FSL naming convention
% 1001_11655_mprage_struc.nii and m1001_11655_mprage_struc_brain.nii

spm('defaults','pet');
spm_jobman('initcfg');

nsegs = 6; %THis is automatically overridden
do_fillholes = 1;  %Fill holes in brainmask
do_clean = 1;  %remove small bits
do_erode = 1;
do_dilate =0;  %Not done for endo t1 k images = 1 for smooth epi mask
medfilt = 5;
connectivity = 6; %or 18 or 26

% Make brain smaller = make non-brain bigger --> non_brain_thresh smaller
% Make brain bigger = non_brain_thresh bigger
brain_thresh = 0.5; %It does not matter what this is as non_brain_thresh controls
non_brain_thresh6seg = 0.1;
non_brain_thresh5seg = 0.3;
non_brain_thresh3seg = 0.3;

se = strel('ball',3,3);  %floor(se+1)/2 with vozel size of [ 1 1 1]

clean_vol = 6000; %mm3 - enough for vitamin marker 

mask_binary = 0;  % If true, make a binary brain mask from segments rather than the T1 masked
t1_ss_name_pre = 'mask_'; %prefix for mask name 'mask_' used for EPI mask

version=spm('Ver');
if nargin<1
    Pt1 = spm_select([1,Inf],'IMAGE','Select T1 images');
end

disp('Preparing data ...');
nsubjects=size(Pt1,1);

for nP=1:nsubjects   

    disp(['Masking image ' num2str(nP) ' of ' num2str(nsubjects) ' subjects ...']);
    
    [inpath,infile,inext]=spm_fileparts(Pt1);
    
    allfiles=dir(inpath);
    
%     %Deal with using bias corrected image
%     if 1
%         infile=infile(2:end);
%     end
    cims = regexpi({allfiles.name},strcat(strcat('^c','[123456]',infile),inext),'match');
    cims = [cims{:}];
    
    %override explicit assignment of number of segments to use
    nsegs=numel(cims);
    Pseg=cell(1,nsegs);
    
    for nseg=1:nsegs       
        Pseg{nseg}=fullfile(inpath,strcat(strcat('c',num2str(nseg),infile),inext));
    end

    Vseg = spm_vol(strvcat(Pseg));
    Vt1 = spm_vol(Pt1);
    
    Yt1 = spm_read_vols(Vt1);
    Yseg = spm_read_vols(Vseg);

    if 2==nsegs
        non_brain=zeros(size(Yseg(:,:,:,1)));        
    elseif 3==nsegs
        non_brain = Yseg(:,:,:,3)  > non_brain_thresh3seg;
    elseif 5==nsegs
        non_brain = Yseg(:,:,:,3) + Yseg(:,:,:,4) + Yseg(:,:,:,5) > non_brain_thresh5seg;
    elseif 6==nsegs
        non_brain = Yseg(:,:,:,3) + Yseg(:,:,:,4) + Yseg(:,:,:,5) + Yseg(:,:,:,6) > non_brain_thresh6seg;
    else
        error(['Is brain segmented? Did not find 2, 3, 5 or 6 tissue segments to match ' Pt1]);
    end
    brain = (Yseg(:,:,:,1) + Yseg(:,:,:,2)) > brain_thresh;    
    brain = brain & ~non_brain;
    brain = single(brain);  %%Stops imfill warning
    
    if do_fillholes
%        for repeat = 1 : 3
            for k=1:Vseg(1).dim(3)
                brain(:,:,k)=imfill(squeeze(brain(:,:,k)),'holes');  % Fill holes in mask
            end
%        end
    end
    
    if do_clean
        voxdims=sqrt(sum(Vt1.mat(1:3,1:3).^2));
        clean_vox = round(clean_vol / prod(voxdims));
        [brain,num]=bwlabeln(brain,connectivity);
        
        stats = regionprops(brain,'Area');
        idx = find([stats.Area] > clean_vox);
        brain = ismember(brain,idx);
%         for j=1:num
%             brain_tmp=brain==j;
%             if sum(brain_tmp(:)) < clean_vox
%                 brain(brain==j)=0;
%             end
%         end
        
%         brain =logical(brain);
    end
    
    if do_erode
            brain=imdilate(single(brain),se);  %Fills holes
            brain=imerode(single(brain),se);  %Fills holes

            brain=brain>0.1;
            brain = single(brain);  %%Stops imfill warning
            for k=1:Vseg(1).dim(3)
                brain(:,:,k)=imfill(brain(:,:,k),'holes');  % Fill holes in mask
            end                                 
    end
 
    if do_dilate
        
        brain=imdilate(single(brain),se);  %Fills holes
        %brain=imerode(single(brain),se);  %Erodes
        
        brain=(brain-min(brain(:)))/(max(brain(:))-min(brain(:)));
        brain=brain>0;
        
        brain = single(brain);  %%Stops imfill warning
        for k=1:Vseg(1).dim(3)
            brain(:,:,k)=imfill(brain(:,:,k),'holes');  % Fill holes in brain
        end            
        %brain=brain>0.5;
        brain = single(brain);  %%Stops imfill warning
        
       for k=1:Vseg(1).dim(3)
            brain(:,:,k)=medfilt2(brain(:,:,k),[medfilt,medfilt]);
            brain(:,:,k)=imfill(brain(:,:,k),'holes');
       end
 
        for k=1:Vseg(1).dim(2)
            brain(:,k,:)=medfilt2(squeeze(brain(:,k,:)),[medfilt,medfilt]);
            brain(:,k,:)=imfill(squeeze(brain(:,k,:)),'holes');
        end
        
        for k=1:Vseg(1).dim(1)
            brain(k,:,:)=medfilt2(squeeze(brain(k,:,:)),[medfilt,medfilt]);
            brain(k,:,:)=imfill(squeeze(brain(k,:,:)),'holes');
        end
%             imdilate(brain, ones(medfilt, medfilt, medfilt));
%          for k=1:Vseg(1).dim(3)
%             brain(:,:,k)=imfill(brain(:,:,k),'holes');  % Fill holes in brain
%          end         
 
        brain=brain~=0;                      
    end     
    
    t1_ss = Yt1 .* brain;

    Vt1_ss = Vt1;
    out_dir=inpath;
    if mask_binary
        t1_ss = brain;    
        t1_ss_name = fullfile(out_dir,strcat(t1_ss_name_pre,infile,inext));
        t1_ss = logical(t1_ss);

        Vt1_ss.dt(1)=spm_type('uint8');
        Vt1_ss.pinfo(1)=1/(2^8-1);
    else
        t1_ss = Yt1 .* brain;
        t1_ss_name = fullfile(out_dir,strcat(infile,'_struc_brain',inext));
        %COPY ORIGINAL files to a specified directory for FSL VBM
%         t1_ss_name = fullfile(out_dir,strcat(aname,'_struc_brain',inext));
        %COPY ORIGINAL files for FSL VBM
%         copy_nii_file_to(Vt1.fname,fullfile(out_dir,strcat(aname,'_struc')));
    end    
    
    Vt1_ss.fname = t1_ss_name;
    Vt1_ss.descrip = ['brain masked by ' nsegs ' segs ' num2str(brain_thresh)];
    
    Vout = spm_create_vol(Vt1_ss);

    for p=1:Vt1_ss.dim(3)
        Vout = spm_write_plane(Vout,t1_ss(:,:,p),p);
    end        
end    


disp('Finished skull stripping data.');



