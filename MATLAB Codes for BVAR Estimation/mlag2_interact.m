function [Ylag] = mlag2_interact(Y,p,W,full)
%MLAG2_interact Create autoregressive lag matrix
%   IF full=true / 1, then
%   This function creates lags of the matrix X(t), in the form:
%       Ylag = [Y(t-1),...,Y(t-p),W(t),W(t-1)Y(t-1),...,W(t-p)Y(t-p)]
%   The interaction is allowed for w and every variable in Y
%
%   IF full=false / 0, then
%   This function creates lags of the matrix X(t), in the form:
%       Ylag = [Y(t-1),...,Y(t-p),W(t),W(t-1)Y_N(t-1),...,W(t-p)Y_N(t-p)]
%   The interaction is only allowed for w and Y_N
%
%   Y is [T x M]. 
%   Ylag with simple interaction is [T x (M*p+w+p*w)]
%   Ylag with full interaction is [T x (M*p+w+M*p*w)]
%
%   Written by Dimitris Korobilis, March 2007
%   Modified by Zizhe Xia, 2018

if full
    [Traw,M]=size(Y);
    [~,w]=size(W);
    Ylag=zeros(Traw,M*p+w+M*p*w);
    for i=1:p
        Ylag(p+1:Traw,(M*(i-1)+1):M*i)=Y(p+1-i:Traw-i,1:M);
    end
    Ylag(1:Traw,M*p+1:M*p+w)=W;
    for k=1:p % iter lag period
        for j=1:w % iter var in w
%             size(Y(p+1-k:Traw-k,1:M))
%             size(W(p+1-k:Traw-k, j:j))
%             W(p+1-k:Traw-k, j:j) .* Y(p+1-k:Traw-k,1:M)
            Ylag(p+1:Traw, M*p+w+1+(j-1)*M+(k-1)*w*M:M*p+w+j*M+(k-1)*w*M) = W(p+1-k:Traw-k, j:j) .* Y(p+1-k:Traw-k,1:M);
        end
    end
%     for i_2=1:p
% %         Ylag(p+1:Traw,N*p+M+1+(i_2-1)*N:N*p+M+1+i_2*N)=w(p+1-i_2:Traw-i_2,:) .* Y(p+1-i_2:Traw-i_2,:);
%         for i_3=1:M
%            Ylag(p+1:Traw,N*p+M+1+(i_2-1)*N*M+(i_3-1)*M:N*p+M+1+i_2*N*M+i_3*M) = w(p+1-i_2:Traw-i_2,i_3) .* Y(p+1-i_2:Traw-i_2,:);
%         end
%     end
else
    [Traw,M]=size(Y);
    [~,w]=size(W);
    Ylag=zeros(Traw,M*p+w+p*w);
    for i=1:p
        Ylag(p+1:Traw,(M*(i-1)+1):M*i)=Y(p+1-i:Traw-i,1:M);
    end
    Ylag(1:Traw,M*p+1:M*p+w)=W;
    for k=1:p
        for j=1:w
            Ylag(p+1:Traw,M*p+w+j+(k-1)*w:M*p+w+j+(k-1)*w)=W(p+1-k:Traw-k,j:j) .* Y(p+1-k:Traw-k,M:M);
        end % M*p+w+1+(j-1)*M+(k-1)*w*M:M*p+w+j*M+k*w*M)
    end
%     for i_2=1:p
%         Ylag(p+1:Traw,N*p+1+i_2:N*p+1+i_2)=w(p+1-i_2:Traw-i_2,1:1) .* Y(p+1-i_2:Traw-i_2,N:N);
%     end
end

% %OR:
% [Traw,N]=size(X);
% Ylag=zeros(Traw,N,p);
% for i_1=1:p
%     Ylag(p+1:Traw,:,i_1)=X(p+1-i_1:Traw-i_1,:);
% end
% Ylag=Ylag(:,:);