classdef wavelet < matlab.mixin.Copyable
    properties
        name
        scale
        transp
        S
    end
    properties ( Access = private )
    end
    methods
        function b = wavelet(name, scale)
            b.name = name;
            b.scale = scale;
            b.transp = 0;
            b.S = [];
        end
        function b = mtimes(obj, a)
            if length(size(a)) == 4 % if multiple frames do for each frame
                if size(a,3) == 1 %2D
                    for t = 1:size(a,4)
                        if obj.transp
                            br = waverec2(real(a(:,:,:,t)), obj.S, obj.name); 
                            bi = waverec2(imag(a(:,:,:,t)), obj.S, obj.name); 
                            b(:,:,:,t)  = br + 1j*bi;
                        else
                            [br, obj.S] = wavedec2(real(a(:,:,:,t)), obj.scale, obj.name);        
                            [bi, obj.S] = wavedec2(imag(a(:,:,:,t)), obj.scale, obj.name);        
                            b(:,:,:,t) = br + 1j*bi;
                        end
                    end
                elseif size(a,3) > 1 %3D
                    for t = 1:size(a,4)
                        if obj.transp
                            br = waverec3(real(a(:,:,:,t)), obj.S, obj.name); 
                            bi = waverec3(imag(a(:,:,:,t)), obj.S, obj.name); 
                            b(:,:,:,t)  = br + 1j*bi;
                        else
                            [br, obj.S] = wavedec3(real(a(:,:,:,t)), obj.scale, obj.name);        
                            [bi, obj.S] = wavedec3(imag(a(:,:,:,t)), obj.scale, obj.name);        
                            b(:,:,:,t) = br + 1j*bi;
                        end
                    end
                end
            else % single frame
                if size(a,3) == 1 %2D
                    if obj.transp
                        br = waverec2(real(a), obj.S, obj.name); 
                        bi = waverec2(imag(a), obj.S, obj.name); 
                        b  = br + 1j*bi;
                    else
                        [br, obj.S] = wavedec2(real(a), obj.scale, obj.name);        
                        [bi, obj.S] = wavedec2(imag(a), obj.scale, obj.name);        
                        b = br + 1j*bi;
                    end
                elseif size(a,3) > 1 %3D
                    if obj.transp
                        br = waverec3(real(a), obj.S, obj.name); 
                        bi = waverec3(imag(a), obj.S, obj.name); 
                        b  = br + 1j*bi;
                    else
                        [br, obj.S] = wavedec3(real(a), obj.scale, obj.name);        
                        [bi, obj.S] = wavedec3(imag(a), obj.scale, obj.name);        
                        b = br + 1j*bi;
                    end
                end
                
            end
            
        end
        function b = ctranspose(obj)
            b   =   copy(obj);
            b.transp = xor(obj.transp, 1); 
        end
    end
end
