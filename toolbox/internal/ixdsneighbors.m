function varargout = ixdsneighbors(A,D,options)

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

% 1    2   4
% 128      8
% 64  32  16

arguments
    A 
    D = []
    options.keepequal = false
end
keepequal  = options.keepequal;
useauxtopo = ~isempty(D);

roff = [-1 0 1 1];
coff = [1 1 1 0];
bpf = uint8([4 8 16 32]);
bpi = uint8([64 128 1 2]);

siz = size(A);

B   = zeros(siz,'uint8');

for r = 1:siz(1)
    for c = 1:siz(2)         

        for n = 1:4
            rn = r+roff(n);
            cn = c+coff(n);
            if rn > siz(1) || rn == 0
                continue
            end
            if cn > siz(2) 
                continue
            end
            P = A(r,c);
            N = A(rn,cn);

            if P > N
                B(r,c)   = B(r,c)+bpf(n);
                continue
            end
            if P < N
                B(rn,cn) = B(rn,cn) + bpi(n);
                continue
            end
            
            if useauxtopo

                AP = (D(r,c));
                AN = (D(rn,cn));
                if isinf(AP)
                    % do nothing
                elseif isinf(AN)
                    % do nothing
                else
                    if AP > AN
                        B(r,c) = B(r,c)+bpf(n);
                        continue
                    end
                    if AP < AN
                        B(rn,cn) = B(rn,cn) + bpi(n);
                        continue
                    end
                end
            end

            if keepequal && N == P
                    B(r,c)   = B(r,c)+bpf(n);
                    B(rn,cn) = B(rn,cn) + bpi(n);
                    continue
            end       
        end
    end
end


% Now get the number of downstream neighbors in each cell
N = zeros(size(B),'uint8');
for r=1:8
    N = N+bitget(B,r,"uint8");
end
N = sum(N(:),'double');

if nargout == 1
    varargout{1} = N;
    return
end


% Now allocate vectors
ix = zeros(N,1);
ixc = zeros(N,1);
clear N

stride = size(B,1);
offset = [-stride-1 -1 stride-1 stride stride+1 1 -stride+1 -stride];
pos1 = 1;

for r = 1:8
    pix = find(bitget(B,r,'uint8'));
    n   = numel(pix);
    pos2 = pos1 + n - 1;

    ix(pos1:pos2)  = pix;
    ixc(pos1:pos2) = pix+offset(r);
    pos1 = pos2+1;
end

varargout{1} = ix;
varargout{2} = ixc;

% Write a function that returns the sum of the neighbor indicators only for
% those neighbors that are equal or lower than the central pixel
