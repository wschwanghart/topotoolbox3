function [IX,ixcoord,res] = getix(parent,him,options)
%GETIX Get linear index of locations in image
%
% Syntax
%
%     IX = getix
%     IX = getix(parentaxes,himage)
%
% Description
%
%     Given axes with an image (e.g. plotted using imagesc or imageschs),
%     this function lets you draw one or more points and returns the linear
%     index of these locations.
%
% Input arguments
%
%     parentaxes   handle to parent axes of the image (default = gca)
%     himage       handle to image. By default, the image is found by
%                  findobj(parent,'type','image')
%
% Output arguments
%
%     IX           Linear index of pixels in image at mapped location.
%
% Example
%
%     DEM = GRIDobj('srtm_bigtujunga30m_utm11.tif');
%     FD  = FLOWobj(DEM);
%     imageschs(DEM);
%     ix = getix;
%     S  = STREAMobj(FD,'chan',ix);
%     plot(S)
%
% See also: imageschs, coord2ind
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 24. June, 2024

arguments
    parent = gca
    him = findobj(gca,'type','image')
    options.color = [0 0.4470 0.7410]
    options.showhints = true
    options.n = inf
end

if isempty(him)
    error('TopoToolbox:image_not_found',"No image found.")
end

if numel(him)>1
    him = him(1);
    warning('TopoToolbox:getix','Multiple images found')
end

M = size(get(him,'CData'), 1);
N = size(get(him,'CData'), 2);

X = get(him,'XData');
Y = get(him,'YData');

if numel(X) == 2 && numel(Y) == 2
    if isa(X,'cell') && isa(Y,'cell')
        X = X{2};
        Y = Y{2};
    else

        X = linspace(X(1),X(2),N);
        Y = linspace(Y(1),Y(2),M);
    end
end

pos = [];
hold on
hfig = gcf;

htitle = title(parent,'','FontSize',8);
parent.TitleHorizontalAlignment = 'left';

% Drawing area = [x,y,w,h]
drawingarea = [min(X) min(Y) max(X)-min(X) max(Y)-min(Y)];

newPos = true;
n      = 1;
while newPos && n<=options.n
    posn = getPos(parent);
    if ~isempty(posn)
        pos = [pos;posn];
    else
        newPos = false;
    end
    n = n+1;
end

if isempty(pos)
    IX = [];
    ixcoord = [];
    res = [];
else
    [IX,ixcoord,res] = coord2ind(X,Y,pos(:,1),pos(:,2));
end


    function pos = getPos(parent)
        if options.showhints
        htitle.String = 'Draw location.';
        end

        hroi = drawcrosshair(parent,'Selected',true,'Color',options.color,...
            'StripeColor','w','LineWidth',0.1,...
            'DrawingArea',drawingarea);
        if isvalid(hfig)
            set(hfig, 'KeyPressFcn', @(src, event) finalizeCrosshair(src, event, hroi));
            if options.showhints
            htitle.String = ['Move location and double click/press enter to finalize.\newline'... 
                             'Close window or press escape to quit.'];
            end
        end
        
        try

            wait(hroi);
        catch
            pos = [];
            return;
        end

        % Is crosshair deleted?
        if ~isvalid(hroi)
            pos = [];
            return
        end

        % If esc is pressed, than there is an empty position
        if isempty(hroi.Position)
            pos = [];
        else
            pos = hroi.Position;
            plot(hroi.Position(1),hroi.Position(2),'ok',...
                'MarkerFaceColor',[1 1 1],'MarkerSize',4)
            delete(hroi)
        end
    end


    function finalizeCrosshair(~, event, hroi)
        % Check if the pressed key is the Return key (Enter)
        if ~isvalid(hroi)
            return
        end

        if strcmp(event.Key, 'return')
            % Get the final position of the crosshair
            hroi.Position;
        elseif strcmp(event.Key,'h')
            if options.showhints
            
                options.showhints = false;
                htitle.String = '';

            else
                options.showhints = true;
            end

        end
    end
end







