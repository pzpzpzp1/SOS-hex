%{
The trilinear map for the hex with given vertices, evaluated at the given points
	V	vertices of hex, 8x3
	p	points to evaluate at, Nx3
returns:
	phi	X,Y,Z values of the map, Nx3
%}
function phi = symbolic_trilinear_map(V, p)
	phi = [(1 - p(:,0)).*(1 - p(:,1)).*V(1,1).*(1 - p(:,2)) + (1 - p(:,0)).*p(:,1).*V(4,1).*(1 - p(:,2)) + p(:,0).*(1 - p(:,1)).*V(5,1).*(1 - p(:,2)) + p(:,0).*p(:,1).*V(8,1).*(1 - p(:,2)) + (1 - p(:,0)).*(1 - p(:,1)).*V(2,1).*p(:,2) + (1 - p(:,0)).*p(:,1).*V(3,1).*p(:,2) + p(:,0).*(1 - p(:,1)).*V(6,1).*p(:,2) + p(:,0).*p(:,1).*V(7,1).*p(:,2), (1 - p(:,0)).*(1 - p(:,1)).*V(1,2).*(1 - p(:,2)) + (1 - p(:,0)).*p(:,1).*V(4,2).*(1 - p(:,2)) + p(:,0).*(1 - p(:,1)).*V(5,2).*(1 - p(:,2)) + p(:,0).*p(:,1).*V(8,2).*(1 - p(:,2)) + (1 - p(:,0)).*(1 - p(:,1)).*V(2,2).*p(:,2) + (1 - p(:,0)).*p(:,1).*V(3,2).*p(:,2) + p(:,0).*(1 - p(:,1)).*V(6,2).*p(:,2) + p(:,0).*p(:,1).*V(7,2).*p(:,2), (1 - p(:,0)).*(1 - p(:,1)).*V(1,3).*(1 - p(:,2)) + (1 - p(:,0)).*p(:,1).*V(4,3).*(1 - p(:,2)) + p(:,0).*(1 - p(:,1)).*V(5,3).*(1 - p(:,2)) + p(:,0).*p(:,1).*V(8,3).*(1 - p(:,2)) + (1 - p(:,0)).*(1 - p(:,1)).*V(2,3).*p(:,2) + (1 - p(:,0)).*p(:,1).*V(3,3).*p(:,2) + p(:,0).*(1 - p(:,1)).*V(6,3).*p(:,2) + p(:,0).*p(:,1).*V(7,3).*p(:,2)]
end