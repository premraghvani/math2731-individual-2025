% seeding for output purposes
rng(1);

% constants
a=4;
b=3;
r=1;
chains = 5;

% function
function [acceptedPhi,acceptedTheta] = convergence(nSamples,burnIn,sigma,firstPhi,firstTheta,a,b,r,chains)
    % init
    acceptedPhi = zeros(nSamples+burnIn,1);
    acceptedTheta = zeros(nSamples+burnIn,1);
    acceptedPhi(1) = firstPhi;
    acceptedTheta(1) = firstTheta;

    % acceptance jacobian
    acceptanceJacobian = @(theta, phi) r .* sqrt( ...
        (sin(theta).^2 .* ( (b + r.*sin(theta)).^2 .* cos(phi).^2 + (a + r.*sin(theta)).^2 .* sin(phi).^2 )) + ...
        (cos(theta).^2 .* ( (b + r.*sin(theta)).*cos(phi).^2 + (a + r.*sin(theta)).*sin(phi).^2 ).^2) ...
    );

    % primary loop
    acceptedElements = 1;
    for i = 2:(nSamples+burnIn)
        proposedPhi = acceptedPhi(i-1) + sigma * randn;
        proposedTheta = acceptedTheta(i-1) + sigma * randn;
    
        % acceptance ratio, choses to accept or reject
        acceptanceRatio = acceptanceJacobian(proposedTheta, proposedPhi) / acceptanceJacobian(acceptedTheta(i-1), acceptedPhi(i-1));
        if rand < acceptanceRatio
            acceptedPhi(i) = proposedPhi;
            acceptedTheta(i) = proposedTheta;
            acceptedElements = acceptedElements+1;
        else
            acceptedPhi(i) = acceptedPhi(i-1);
            acceptedTheta(i) = acceptedTheta(i-1);
        end

        % modulo, eliminate the periodicity
        acceptedPhi(i) = mod(acceptedPhi(i),2*pi);
        acceptedTheta(i) = mod(acceptedTheta(i),2*pi*chains);
    end
    disp(acceptedElements/(nSamples+burnIn))

    % rids of burn in
    acceptedPhi = acceptedPhi(burnIn+1:end);
    acceptedTheta = acceptedTheta(burnIn+1:end);
end

% plot function
function plotter(nSamples,burnIn,sigma,firstPhi,firstTheta,a,b,r,chains,animate)
    [acceptedPhi,acceptedTheta] = convergence(nSamples,burnIn,sigma,firstPhi,firstTheta,a,b,r,chains);

    % label scales, from main
    set(0,"DefaultAxesFontSize",16);
    set(0,"DefaultTextFontSize",18);
    set(0,"DefaultAxesTitleFontSizeMultiplier",1.2);
    set(0,"DefaultAxesLabelFontSizeMultiplier",1.2);
    set(0,"DefaultLineLineWidth",1.5);
    set(0,"DefaultLegendFontSize",12);

    % figure setup
    figure;
    tiledlayout(1,3,"TileSpacing","compact","Padding","compact");

    % plot 1+2: phi and theta
    nbins = 50;
    edges = linspace(0, 2*pi, nbins+1);
    edgesTheta = linspace(0,2*chains*pi,5*nbins+1);
    
    [counts, phiEdges, thetaEdges] = histcounts2(acceptedPhi, acceptedTheta, edges, edgesTheta);
    
    % prior codes used, with copilot auto-do
    dphi = phiEdges(2) - phiEdges(1);
    dtheta = thetaEdges(2) - thetaEdges(1);
    pdf2D = counts / (sum(counts, "all") * dphi * dtheta);
    
    nexttile;
    imagesc(phiEdges, thetaEdges, pdf2D');    
    axis square
    xlabel('\phi');
    ylabel('\theta');
    
    colorbar
    caxis([0 0.1]); % gpt after asking, used in other parts!
    colormap(turbo);
    title(sprintf("$\\phi$ vs $\\theta$ Density Heatmap\n$n = %d,\\; \\sigma = %.2f$", ...
                  nSamples, sigma), 'Interpreter', 'latex');
    
    
    % gets torus density data
    [ZposNormal,ZnegNormal,ZposRotated,ZnegRotated,NposNormal,NnegNormal,NposRotated,NnegRotated,bins,binsX,minZ,maxZ] = torusDensity(acceptedTheta,acceptedPhi,nSamples,a,b,r,chains);

    % plot 3: to emphasise on heights
    nexttile;
    hold on;

    % gridding
    dimMax = max(a,b);
    xe = linspace(-(dimMax+r)*chains, (dimMax+r)*chains, binsX);
    ye = linspace(-(dimMax+r)*1.2, (dimMax+r)*1.2, bins);
    xc = xe(1:end-1) + 0.5*(dimMax+r)/binsX;
    yc = ye(1:end-1) + 0.5*(dimMax+r)/bins;
    [xgrid, ygrid] = meshgrid(xc, yc);
    h3=gca;
    surf(xgrid, ygrid, ZposNormal','EdgeColor','none'); % top unrotated
    surf(xgrid, ygrid, ZnegNormal','EdgeColor','none'); % bottom unrotated
    surf(xgrid, ygrid, ZposRotated','EdgeColor','none'); % top rotated
    surf(xgrid, ygrid, ZnegRotated','EdgeColor','none'); % bottom rotated

    axis equal tight;
    colormap(h3, turbo); colorbar; % gpt - bcause colormap turbo for example is global
    caxis([minZ maxZ]);
    xlabel("X"); ylabel("Y"); zlabel("Z");
    title(sprintf(["Chained Tori\n3D Mapping of Densities"]),"Interpreter", "latex");
    view(45,35);
    camlight headlight; lighting gouraud;

    % plot 4: x,y
    Ntotal = NposNormal + NnegNormal + NposRotated + NnegRotated;
    Ntotal(Ntotal==0) = NaN; 
    nexttile;
    h4 = gca;

    minX = -(a+r); maxX = (a+r)*(chains*2-1); minY = -(b+r); maxY = (b+r);
    % uses old analytic results for a=3,b=2,r=1, but it works fine (seemingly)

    xCoords = linspace(minY,maxY,bins-1);
    yCoords = linspace(minX,maxX,binsX-1); % swapped these around and it works?! - no idea why, but it paints a pretty picture

    hImg = imagesc(xCoords, yCoords, Ntotal); % directly inputted numerically from min, max, interX and Y
    axis equal tight off;
    xlabel("Y"); ylabel("X");
    title(sprintf(["Torus\nDensity Map for X,Y\nAxes Transposed"]),"Interpreter", "latex");
    colormap(h4, hot); colorbar;
    caxis([0 nSamples*0.0003]);

    % chatgpt suggestion to make tiles transparent if NaN - also made hImg
    % up for me
    hImg.AlphaData = ~isnan(Ntotal); 

    % animate procedure
    if animate == true
        animation(acceptedPhi,acceptedTheta,nSamples,sigma,a,b,r,chains)
    end
end

% torus density
function [ZposNormal,ZnegNormal,ZposRotated,ZnegRotated,NposNormal,NnegNormal,NposRotated,NnegRotated,bins,binsX,minZ,maxZ] = torusDensity(acceptedTheta,acceptedPhi,nSamples,a,b,r,chains)
    % torus cartesian formation
    R=max(a,b);
    centerSeparation = 2*(R-r);
    
    chainId = floor(acceptedTheta / (2*pi));
    rotation = mod(chainId,2); % 0 = normal, 1 = side

    x = (a+r*sin(acceptedTheta)).*cos(acceptedPhi) + centerSeparation.*chainId;

    % vectorised variant, as if it is rotated, y and z are flipped
    y = (b+r*sin(acceptedTheta)).*sin(acceptedPhi).*abs(rotation-1) + r*cos(acceptedTheta).*rotation;
    z = r*cos(acceptedTheta).*abs(rotation-1) + -1*(b+r*sin(acceptedTheta)).*sin(acceptedPhi).*rotation;


    % inspiration https://uk.mathworks.com/help/matlab/creating_plots/types-of-matlab-plots.html

    % space and bins for x,y axis
    bins = 200;
    binsX = 200*chains;
    
    % logic: we find, in each bin (area), the average of positive z, and average of negative z to shape the torus.
    % can also be used for plot 4 in densities for top view.

    minX = round(min(x),0);
    maxX = round(max(x),0);
    minY = round(min(y),0);
    maxY = round(max(y),0);
    minZ = round(min(z),0);
    maxZ = round(max(z),0);

    binXDivider = (maxX - minX) / binsX;
    binYDivider = (maxY - minY) / bins;

    % number and total Z
    ZposNormal = nan(binsX-1,bins-1);  NposNormal = zeros(binsX-1,bins-1);
    ZnegNormal = nan(binsX-1,bins-1);  NnegNormal = zeros(binsX-1,bins-1);
    ZposRotated = nan(binsX-1,bins-1);  NposRotated = zeros(binsX-1,bins-1);
    ZnegRotated = nan(binsX-1,bins-1);  NnegRotated = zeros(binsX-1,bins-1);
    

    % looping through all items
    for k = 1:nSamples

        % bin indices for this sample
        boxX = floor((x(k) - minX) / binXDivider) + 1;
        boxY = floor((y(k) - minY) / binYDivider) + 1;
        % clip indices
        boxX = min(max(boxX, 1), binsX-1);
        boxY = min(max(boxY, 1), bins-1);

        % positive Z unrotated
        if z(k) >= 0 && rotation(k) == 0
            if isnan(ZposNormal(boxX, boxY))
                ZposNormal(boxX, boxY) = z(k);
            else
                ZposNormal(boxX, boxY) = ZposNormal(boxX, boxY) + z(k);
            end
            NposNormal(boxX, boxY) = NposNormal(boxX, boxY) + 1;

        % positive Z rotated
        elseif z(k) >= 0 && rotation(k) == 1 
            if isnan(ZposRotated(boxX, boxY))
                ZposRotated(boxX, boxY) = z(k);
            else
                ZposRotated(boxX, boxY) = ZposRotated(boxX, boxY) + z(k);
            end
            NposRotated(boxX, boxY) = NposRotated(boxX, boxY) + 1;
        
        % negative Z unrotated
        elseif z(k) < 0 && rotation(k) == 0
            if isnan(ZnegNormal(boxX, boxY))
                ZnegNormal(boxX, boxY) = z(k);
            else
                ZnegNormal(boxX, boxY) = ZnegNormal(boxX, boxY) + z(k);
            end
            NnegNormal(boxX, boxY) = NnegNormal(boxX, boxY) + 1;

        % negative Z rotated
        else
            if isnan(ZnegRotated(boxX, boxY))
                ZnegRotated(boxX, boxY) = z(k);
            else
                ZnegRotated(boxX, boxY) = ZnegRotated(boxX, boxY) + z(k);
            end
            NnegRotated(boxX, boxY) = NnegRotated(boxX, boxY) + 1;
        end
    end

    % finds average density per box
    ZposNormal = ZposNormal ./ NposNormal;
    ZnegNormal = ZnegNormal ./ NnegNormal;
    ZposRotated = ZposRotated ./ NposRotated;
    ZnegRotated = ZnegRotated ./ NnegRotated;
    ZposNormal(NposNormal==0) = nan;
    ZnegNormal(NnegNormal==0) = nan;
    ZposRotated(NposRotated==0) = nan;
    ZnegRotated(NnegRotated==0) = nan;
end

% animation func
function animation(acceptedPhi,acceptedTheta,nSamples,sigma,a,b,r,chains)
    nFrames = 450;
    fig = figure('Position', [0 50 1920 1080]); % new tab
    tiledlayout(fig, 1, 3, "TileSpacing","compact","Padding","compact");

    % starts video
    v = VideoWriter('torus_animation_elliptic_chained.mp4','MPEG-4');
    v.FrameRate = 25;
    open(v);

    % theta phi map
    nexttile;
    nbins = 50;
    edges = linspace(0, 2*pi, nbins+1);
    edgesTheta = linspace(0,2*chains*pi,5*nbins+1);
    [counts, phiEdges, thetaEdges] = histcounts2(acceptedPhi(1:1),acceptedTheta(1:1),edges,edgesTheta);

    dphi = phiEdges(2)-phiEdges(1);
    dtheta = thetaEdges(2)-thetaEdges(1);
    pdf2D = counts / (sum(counts,"all") * dphi * dtheta);

    hHeat = imagesc(phiEdges,thetaEdges,pdf2D');
    axis square
    xlabel("\phi"); ylabel("\theta");
    colorbar

    caxis([0 0.1]);
    xlim([0 2*pi]);
    ylim([0 2*chains*pi]);

    colormap(turbo);

    % theta (removed – now absorbed into heatmap)
    % (no second histogram tile)

    % key torus values
    [ZposNormal,~,~,~,NposNormal,~,~,~,bins,binsX,~,~]= torusDensity(acceptedTheta(1:1), acceptedPhi(1:1), 1,a,b,r,chains );

    % torus 3d
    nexttile;
    dimMax = max(a,b);
    xe = linspace(-(dimMax+r)*chains, (dimMax+r)*chains, binsX);
    ye = linspace(-(dimMax+r)*1.2, (dimMax+r)*1.2, bins);
    xc = xe(1:end-1) + 0.5*(dimMax+r)/binsX;
    yc = ye(1:end-1) + 0.5*(dimMax+r)/bins;
    [xgrid3D, ygrid3D] = meshgrid(xc, yc);

    h3 = gca;
    hold on;
    ZposPlot = surf(xgrid3D, ygrid3D, nan(size(ZposNormal')),'EdgeColor','none'); % top
    ZnegPlot = surf(xgrid3D, ygrid3D, nan(size(ZposNormal')),'EdgeColor','none'); % bottom
    ZposRotatedPlot = surf(xgrid3D, ygrid3D, nan(size(ZposNormal')),'EdgeColor','none'); % top rotate
    ZnegRotatedPlot = surf(xgrid3D, ygrid3D, nan(size(ZposNormal')),'EdgeColor','none'); % bottom rotate

    axis equal tight off;
    colormap(h3, turbo); colorbar;
    
    caxis([-(b+r) (b+r)]); % numeric known
    xlim([xe(1) xe(end)])
    ylim([ye(1) ye(end)])
    zlim([-(b+r) (b+r)])

    view(45,35); camlight headlight; lighting gouraud;

    % birds eye
    nexttile;
    h4 = gca;

    minX = -(a+r); maxX = (a+r)*(chains*2-1); minY = -(b+r); maxY = (b+r);
    xCoords = linspace(minY,maxY,bins-1);
    yCoords = linspace(minX,maxX,binsX-1); % swapped these around and it works?! - no idea why, but it paints a pretty picture

    hImg = imagesc(xCoords, yCoords, NposNormal);
    hImg.AlphaData = ~isnan(NposNormal);
    axis equal tight off;

    caxis([0 nSamples*0.0003]);
    %xlim([-(R+r+0.05) (R+r+0.05)])
    %ylim([-(R+r+0.05) (R+r+0.05)])

    colormap(h4, hot); colorbar;

    % summary msg
    msg = sprintf('sigma = %.3f\nn = %d', sigma, nSamples);
    
    note = annotation(fig, 'textbox', ...
    [0.01, 0.01, 0.25, 0.1], ...
    "String", msg, ...
    "FontSize", 14, ...
    "HorizontalAlignment", "left", ...
    "VerticalAlignment", "bottom", ...
    "BackgroundColor", "none", ...
    "EdgeColor", "none", ...
    "Color", "white");  


    % frame update
    frameIdx = unique( round( logspace(0, log10(nSamples), nFrames) ) );
    for i = frameIdx
        [counts, ~, ~] = histcounts2( ...
            acceptedPhi(1:i), acceptedTheta(1:i), edges, edges);
        pdf2D = counts / (sum(counts,"all") * dphi * dtheta);
        set(hHeat,"CData",pdf2D');

        % recompute torus density
        [ZposNormal,ZnegNormal,ZposRotated,ZnegRotated,NposNormal,NnegNormal,NposRotated,NnegRotated,~,~,~,~] = torusDensity(acceptedTheta(1:i), acceptedPhi(1:i), i,a,b,r,chains );
        set(ZposPlot, "ZData",ZposNormal');
        set(ZnegPlot,"ZData",ZnegNormal');
        set(ZposRotatedPlot, "ZData",ZposRotated');
        set(ZnegRotatedPlot,"ZData",ZnegRotated');

        Ntotal = NposNormal + NnegNormal + NposRotated + NnegRotated;
        Ntotal(Ntotal==0) = NaN; 
        hImg.CData = Ntotal;
        hImg.AlphaData = ~isnan(Ntotal);

        note.String = sprintf('sigma = %.3f\nn = %d / %d', sigma, i, nSamples);

        drawnow limitrate   
        frame = getframe(fig);
        writeVideo(v, frame); 
    end

    close(v);
end


% primary seq
plotter(1e8,1e4,0.25,pi,pi,a,b,r,chains,true)
%plotter(1e5,1e4,0.25,pi,pi,a,b,r,chains,false)