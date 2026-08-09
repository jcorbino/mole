function tgvPlotFrame(t, s, step, steps, X, Y, m, n, a, b, c, d)
% Callback for taylor_green_vortex.m: draws the scalar field every few steps.
%
% Lives in its own file because a local function cannot satisfy both languages
% at once -- MATLAB requires script-local functions at the end of the file,
% Octave requires them before the point of use.

    if mod(step, 5) ~= 0 && step ~= steps
        return
    end

    pcolor(X', Y', reshape(s, m+2, n+2)')
    shading interp
    axis equal
    axis([a b c d])
    clim([0 max(s)])   % rescaled every frame; diffusion drops the peak ~4x
    title(sprintf('t = %1.2f s', t))
    xlabel('x')
    ylabel('y')
    colormap jet
    colorbar
    set(gcf, 'color', 'w')
    drawnow
end
