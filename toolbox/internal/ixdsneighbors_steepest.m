function varargout = ixdsneighbors_steepest(A,D,cs)

%IXDSNEIGHBORS Return list of edges connecting pixels to downstream neighbors
%
% Syntax
%
%     [ix,ixc] = ixdsneighbors(A)
%     N = ixdsneighbors(A)
%     ... = ixdsneighbors(A,keepequal = false)
%     ... = ixdsneighbors(A,D)
%
% Description
%
%     The function returns the edge list [ix,ixc] that contains the linear
%     indices that connect the pixels ix with their downstream neighbors
%     ixc. 
%
%     When called with one output argument, the function returns the number
%     of edges.
%
%     The keepequal option determines whether pixel links should be
%     established for pixels with equal height. In flat areas, the option
%     leads to an undirected network.
%
%     To calculate a directed network in flat areas, provide another grid
%     with auxiliary topography D. Elements in D have value -inf unless
%     they are located in flat areas. Here, an auxiliary topography will be
%     used. 
%
% See also: FLOWobj
%
% Author: Wolfgang Schwanghart (schwangh[at]uni-potsdam.de)
% Date: 6. September, 2024

% 1    2   4
% 128      8
% 64  32  16

% 1   2   3
% 8       4
% 7   6   5

arguments
    A 
    D 
    cs
end
useauxtopo = ~isempty(D);

roff = [0 1 -1 1];
coff = [1 0 1  1];

bpf = uint8([4 6 3 5]);
bpi = uint8([8 2 7 1]);
cs2 = sqrt(2*cs);
dx  = [cs cs cs2 cs2];

siz = size(A);
G   = zeros(siz,'single');
B   = zeros(siz,'uint8');

for r = 1:siz(1)
    for c = 1:siz(2)   
        P = A(r,c);
        if isnan(P)
            continue
        end

        for n = 1:4
            rn = r+roff(n);
            cn = c+coff(n);
            if rn > siz(1) || rn == 0
                continue
            end
            if cn > siz(2) 
                continue
            end

            N = A(rn,cn);
            if isnan(N)
                continue
            end
               
            % Calculate slope
            s = (P-N)/dx(n);

            if s > 0
                if s > G(r,c)
                    G(r,c) = s;
                    B(r,c) = bpf(n);
                end
                continue
            end
            if s < 0
                s = -s;
                if s > G(rn,cn)
                    G(rn,cn) = s;
                    B(rn,cn) = bpi(n);
                end
                continue
            end
            
            if useauxtopo

                AP = (D(r,c));
                AN = (D(rn,cn));
                if isinf(AP)
                    continue
                elseif isinf(AN)
                    continue
                else
                    s = (AP-AN)/dx(n);

                    if s > 0
                        if s > G(r,c)
                            G(r,c) = s;
                            B(r,c) = bpf(n);
                        end
                        continue
                    end
                    if s < 0
                        s = -s;
                        if s > G(rn,cn)
                            G(rn,cn) = s;
                            B(rn,cn) = bpi(n);
                        end
                        continue
                    end
                    
                end
            end     
        end
    end
end

stride = size(B,1);
offset = [-stride-1 -1 stride-1 stride stride+1 1 -stride+1 -stride];
ix  = find(B);
ixc = ix + double(offset(B(ix)))';
varargout{1} = uint32(ix);
varargout{2} = uint32(ixc);

