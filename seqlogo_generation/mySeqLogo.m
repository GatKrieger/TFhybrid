function [] = mySeqLogo(pfm, varargin)
% from Sagie Brodsky
p = inputParser();
p.addParameter('subplotDirections', []);
p.parse(varargin{:});
spDir = p.Results.subplotDirections;

[wm,~] = seqlogo(pfm,'DisplayLogo',false);
letters = wm{1};
clr = [0, 153, 51; 0, 45, 179; 204, 153, 0 ; 153, 0, 0]./255; % corresponding colors
%     clr = [0, 77, 0 ; 64, 0, 128;  204, 153, 0 ; 153, 0, 0]./255;
for t = 1:numel(letters)
    if exist([letters(t) '.png'], 'file')==2
        continue;
    end
    hf = figure('position',[100 100 800 810],'color','w');
    ha = axes('parent',hf, 'visible','off','position',[0 0 1 1]);
    ht = text(400,450,letters(t),'color',clr(t,:),'units','pixels',...
        'fontsize',800,'fontweight','norm',...
        'vertical','mid','horizontal','center');
    
    F = getframe(hf); % rasterize the letter
    img = F.cdata;
    m = any(img < 255,3); % convert to binary image
    m(any(m,2),any(m,1))=1; % mask to cut white borders
    imwrite(reshape(img(repmat(m,[1 1 3])),[sum(any(m,2)) sum(any(m,1)) 3]),...
        [letters(t) '.png'])
    close(hf)
end

xlabels = cellstr(num2str([1:size(wm{2},2)]'));
%letters = wm{1};
letters = ['A','C','G','T'];
wmat=wm{2}; % weight matrix from seqlogo
[nletters  npos] = size(wmat);
wmat(wmat<0) = 0; % cut negative values

% prepare the figure
if isempty(spDir)
    f = figure;
else
    subplot(spDir(1), spDir(2), spDir(3))
end
% hAx = axes('parent',gcf,'visible','on');
% hAx = axes;
% set(hAx,'XLim',[0.5 npos+0.5],'XTick',1:npos,'XTickLabel',xlabels)
set(gca,'XLim',[0.5 npos+0.5],'XTick',1:npos,'XTickLabel',xlabels)
ymax = ceil(max(sum(wmat)));
ylim([0 ymax])
%axpos = get(hAx,'Position');
axpos = get(gca,'Position');
step = axpos(3)/npos;

% place images of letters
for i=1:npos
    [wms idx] = sort(wmat(:,i)); % largest on the top
    let_show = letters(idx);
    ybot = axpos(2);
    for s=1:nletters
        if wms(s)==0, continue, end;
        axes('position',[axpos(1) ybot step wms(s)/ymax*axpos(4)])
        ybot = ybot + wms(s)/ymax*axpos(4);
        img = imread([let_show(s) '.png']);
        image(img)
        set(gca,'visible','off')
    end
    axpos(1)=axpos(1)+step;
end
set(gcf, 'color', 'w');
end

  