% seeding for output purposes
rng(1);

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
function plotter(nSamples,burnIn,sigma,firstPhi,firstTheta)
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

    % plot 1: phi
    nexttile;
    h=histogram(acceptedPhi, 50,"Normalization", "pdf", "DisplayName", "Generated Sample","BinLimits",[0,2*pi]);
    xlim([0,2*pi]);
    ylim([0,1.5/(2*pi)]);
    axis square 
    xlabel("\phi")
    ylabel("Frequency")
    hold on;
    xgrid = linspace(0,2*pi,1000);
    y = ones(size(xgrid)) * (1 / (2*pi));
    plot(xgrid, y, 'LineWidth', 2, 'DisplayName', 'Constant $\frac{1}{2\pi}$');
    legend show
    legend("Location","southwest");
    legend("Interpreter","latex");
    title(sprintf(["$\\phi$\nBurned in $\\phi_0 = %.4f$\n$n = %d,\\; \\sigma = %.2f$"],acceptedPhi(1), nSamples, sigma),"Interpreter", "latex");


    % plot 2: theta
    nexttile;
    h=histogram(acceptedTheta, 50,"Normalization", "pdf", "DisplayName", "Generated Sample","BinLimits",[0,2*pi]);
    xlim([0,2*pi]);
    ylim([0,1.5/(2*pi)]);
    axis square 
    xlabel("\theta")
    ylabel("Frequency")
    hold on;
    xgrid = linspace(0,2*pi,1000);
    y = ones(size(xgrid)) * (1 / (2*pi));
    plot(xgrid, y, 'LineWidth', 2, 'DisplayName', 'Constant $\frac{1}{2\pi}$');
    legend show
    legend("Location","southwest");
    legend("Interpreter","latex");
    title(sprintf(["$\\theta$\nBurned in $\\theta_0 = %.4f$\n$n = %d,\\; \\sigma = %.2f$"],acceptedTheta(1), nSamples, sigma),"Interpreter", "latex");

    % torus cartesian formation
    R=2;
    r=1;
    x = (R+r*sin(acceptedTheta)).*cos(acceptedPhi);
    y = (R+r*sin(acceptedTheta)).*sin(acceptedPhi);
    z = r*cos(acceptedTheta);

    % inspiration https://uk.mathworks.com/help/matlab/creating_plots/types-of-matlab-plots.html

    % space and bins for x,y axis
    bins = 200;
    xe = linspace(-R*1.5, R*1.5, bins);
    ye = linspace(-R*1.5, R*1.5, bins);
    
    % logic: we find, in each bin (area), the average of positive z, and average of negative z to shape the torus.
    % can also be used for plot 4 in densities for top view.

    % number and total Z
    Zpos = nan(bins-1,bins-1);  Npos = zeros(bins-1,bins-1);
    Zneg = nan(bins-1,bins-1);  Nneg = zeros(bins-1,bins-1);
    

    % looping through all items
    for k = 1:nSamples

        % bin indices for this sample
        boxX = floor((x(k) + R + r) * bins / (2*(R+r)));
        boxY = floor((y(k) + R + r) * bins / (2*(R+r)));
        % clip indices
        boxX = min(max(boxX, 1), bins-1);
        boxY = min(max(boxY, 1), bins-1);

        % positive Z
        if z(k) >= 0
            if isnan(Zpos(boxX, boxY))
                Zpos(boxX, boxY) = z(k);
            else
                Zpos(boxX, boxY) = Zpos(boxX, boxY) + z(k);
            end
            Npos(boxX, boxY) = Npos(boxX, boxY) + 1;
        else % negative Z
            if isnan(Zneg(boxX, boxY))
                Zneg(boxX, boxY) = z(k);
            else
                Zneg(boxX, boxY) = Zneg(boxX, boxY) + z(k);
            end
            Nneg(boxX, boxY) = Nneg(boxX, boxY) + 1;
        end
    end

    % finds average density per box
    Zpos = Zpos ./ Npos;
    Zneg = Zneg ./ Nneg;
    Zpos(Npos==0) = nan;
    Zneg(Nneg==0) = nan;
    
    
    % plot 3: to emphasise on heights
    nexttile;
    hold on;

    % gridding
    xc = xe(1:end-1) + 0.5*(R+r)/bins;
    yc = ye(1:end-1) + 0.5*(R+r)/bins;
    [xgrid, ygrid] = meshgrid(xc, yc);
    h3=gca;
    surf(xgrid, ygrid, Zpos','EdgeColor','none'); % top
    surf(xgrid, ygrid, Zneg','EdgeColor','none'); % bottom

    axis equal tight;
    colormap(h3, turbo); colorbar; % chatgpt - bcause colormap turbo for example is global
    xlabel("X"); ylabel("Y"); zlabel("Z");
    title(sprintf(["Torus\n3D Mapping of Densities"]),"Interpreter", "latex");
    view(45,35);
    camlight headlight; lighting gouraud;

    % plot 4: x,y
    Ntotal = Npos + Nneg;
    Ntotal(Ntotal==0) = NaN; 
    nexttile;
    h4 = gca;
    hImg = imagesc(-(R+r)+(R+r)/bins:2*(R+r)/bins:(R+r), -(R+r)+(R+r)/bins:2*(R+r)/bins:(R+r), Ntotal);
    axis equal tight;
    xlabel("X"); ylabel("Y");
    title(sprintf(["Torus\nDensity Map for X,Y"]),"Interpreter", "latex");
    colormap(h4, hot); colorbar;
    
    % chatgpt suggestion to make tiles transparent if NaN - also made hImg
    % up for me
    hImg.AlphaData = ~isnan(Ntotal); 
end

plotter(1e7,1e4,0.25,pi,pi)
plotter(1e5,1e4,0.25,pi,pi)
plotter(1e7,1e4,0.005,pi,pi)