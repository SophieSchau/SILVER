function [] = makegif_fast(filename)


      % Capture the plot as an image 
      frame = getframe(gcf); 
      im = frame2im(frame); 
      [imind,cm] = rgb2ind(im,256); 
      % Write to the GIF File 
      if ~exist(filename,'file') 
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime',0.1); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime',0.1); 
      end 
end

