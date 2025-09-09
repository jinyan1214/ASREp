function [ux,uy,uz] = u_3D_comas_v2(Vltp, Kx, Ky, x, y, z, ...
    delta, z0, ys, yf, d, groundNodes)
    
    Smax = (Vltp/100 * pi * d^2)/sqrt(2*pi)/Kx/(z0+min(min(z)))/4;
    ux = zeros(size(x,2),1);
    uy = zeros(size(x,2),1);
    uz = zeros(size(x,2),1);
    y0 = -1 * norminv(delta) * Ky * z0;
    
    x = x(groundNodes)';
    y = y(groundNodes)';
    z = z(groundNodes)';
    uz(groundNodes) = -1 * Smax * exp(-x.^2/2/Kx^2./(z0-z).^2) ...
        .* (normcdf((y-ys-y0)/Ky./(z0-z)) - ...
        normcdf((y-yf)/Ky./(z0-z)));
    
    ux(groundNodes) = x./(z0-z).*uz(groundNodes);
    uy(groundNodes) = Vltp/100*d^2/8./(z0-z).*...
        (exp((-(y-ys-y0).^2 - x.^2)/2/Ky^2./(z0-z).^2) - ...
        exp((-(y-yf).^2 - x.^2)/2/Ky^2./(z0-z).^2));
    
    
end