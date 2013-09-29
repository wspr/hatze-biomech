function patch_rotate(x,y,z,O,R,varargin)

[f,v,c] = surf2patch(x,y,z);

rv = repmat(O,[1 length(v)]) + R*transpose(v);

patch(struct('Faces',f,'Vertices',transpose(rv)),varargin{:});

end