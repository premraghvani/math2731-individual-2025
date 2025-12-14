% WARNING:
% Unlike the main, which considers theta up to 10 pi, this only considers
% up to 2pi, and separately determines which torus (0 to 4). Although WLLN
% means it should be fine, there will be some problems!

% seeding for output purposes
rng(4);

% function
function [acceptedPhi,acceptedTheta] = convergence(nSamples,burnIn,sigma,firstPhi,firstTheta)
    % init
    acceptedPhi = zeros(nSamples+burnIn,1);
    acceptedTheta = zeros(nSamples+burnIn,1);
    acceptedPhi(1) = firstPhi;
    acceptedTheta(1) = firstTheta;

    % primary loop
    for i = 2:(nSamples+burnIn)
        proposedPhi = acceptedPhi(i-1) + sigma * randn;
        proposedTheta = acceptedTheta(i-1) + sigma * randn;
        acceptedPhi(i) = mod(proposedPhi,2*pi);
        acceptedTheta(i) = mod(proposedTheta,2*pi);
    end

    % rids of burn in
    acceptedPhi = acceptedPhi(burnIn+1:end);
    acceptedTheta = acceptedTheta(burnIn+1:end);
end

% plot function
function plotter(nSamples,burnIn,sigma,firstPhi,firstTheta,animate)
    [acceptedPhi,acceptedTheta] = convergence(nSamples,burnIn,sigma,firstPhi,firstTheta);

    % label scales, from main
    set(0,"DefaultAxesFontSize",16);
    set(0,"DefaultTextFontSize",18);
    set(0,"DefaultAxesTitleFontSizeMultiplier",1.2);
    set(0,"DefaultAxesLabelFontSizeMultiplier",1.2);
    set(0,"DefaultLineLineWidth",1.5);
    set(0,"DefaultLegendFontSize",12);

    % figure setup
    figure;
    tiledlayout(1,4,"TileSpacing","compact","Padding","compact");

    % plot 1: phi and theta
    nbins = 50;
    edges = linspace(0, 2*pi, nbins+1);
    edgesTheta = linspace(0,2*pi,nbins+1);
    
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
    [ZposNormal,ZnegNormal,ZposRotated,ZnegRotated,NposNormal,NnegNormal,NposRotated,NnegRotated,R,r,bins,binsX,minZ,maxZ,chainId] = torusDensity(acceptedTheta,acceptedPhi,nSamples);

    % plot 2: histogram of the chains
    nexttile;
    counts = histcounts(chainId, -0.5:1:4.5);
    bar(0:4, counts);
    axis square
    hold on;
    yline(nSamples/5,"LineStyle","--","LineWidth",2,"Color","r");
    xlabel("Chain ID (0=Left, 4=Right)");
    ylabel("Count");
    
    title(sprintf(["Bias Between Chains"]),"Interpreter", "latex");
    xticks(0:4);

    % plot 3: to emphasise on heights
    nexttile;
    hold on;

    % gridding
    xe = linspace(-(R+r)*5, (R+r)*5, binsX);
    ye = linspace(-(R+r)*1.2, (R+r)*1.2, bins);
    xc = xe(1:end-1) + 0.5*(R+r)/binsX;
    yc = ye(1:end-1) + 0.5*(R+r)/bins;
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

    minX = -4; maxX = 20; minY = -3; maxY = 3;
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
end

% torus density
function [ZposNormal,ZnegNormal,ZposRotated,ZnegRotated,NposNormal,NnegNormal,NposRotated,NnegRotated,R,r,bins,binsX,minZ,maxZ,chainId] = torusDensity(acceptedTheta,acceptedPhi,nSamples)
    % torus cartesian formation
    a=4;
    b=3;
    r=1;
    R=max(a,b);
    centerSeparation = 2*(R-r);
    
    % chain ID is a random vector
    chainId = randi([0 4], length(acceptedPhi), 1);

    rotation = mod(chainId,2); % 0 = normal, 1 = side

    x = (a+r*sin(acceptedTheta)).*cos(acceptedPhi) + centerSeparation.*chainId;

    % vectorised variant, as if it is rotated, y and z are flipped
    y = (b+r*sin(acceptedTheta)).*sin(acceptedPhi).*abs(rotation-1) + r*cos(acceptedTheta).*rotation;
    z = r*cos(acceptedTheta).*abs(rotation-1) + -1*(b+r*sin(acceptedTheta)).*sin(acceptedPhi).*rotation;


    % inspiration https://uk.mathworks.com/help/matlab/creating_plots/types-of-matlab-plots.html

    % space and bins for x,y axis
    bins = 200;
    binsX = 800;
    
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

% primary seq
plotter(1e5,1e4,0.25,pi,pi,false)