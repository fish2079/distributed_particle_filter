function beta = multiVarConv(alpha)
%This function calculates "convolution of coefficients of multi-variate 
%polynomials". See Eq. 8 of our Asilomar 2010 paper. This function however
%uses the double sum notation used in pre-final versions of the paper.
%
%alpha - vector of coefficients of polynomial approximation of measurement
%function
%
%beta - 1 x numJLFsums-1 vector of coefficients that need to be summed up 
%using consensus to obtain the coefficients of the joint likelihood.

global parset;

persistent convIndices; %Cell array: There is 1 cell for each entry of beta. This cell contains indices of alphas that need to be multiplied together to obtain the corresponding beta.

numJLFsums = size(parset.powersJLF,1);

if isempty(convIndices) %we execute this only once (during the first call of this function)

    numAlphas = length(alpha);
    
%     for rr = 2:numJLFsums %we do not calculate the beta corresponding to the constant term in the JLF polynomial...
%         r_vect = parset.powersJLF(rr,:); %vector of exponents in the monomial
%         index = 1;
%         for ii = 1:numAlphas
%         i_vect = parset.powersMeasFunc(ii,:); %vector of exponents in the monomial
%             for jj = 1:numAlphas
%                 j_vect = parset.powersMeasFunc(jj,:); %vector of exponents in the monomial
%                 if ((i_vect + j_vect) == r_vect) %check if sum i + j == r...
%                     convIndices{rr}.alpha_i(index) = ii;
%                     convIndices{rr}.alpha_j(index) = jj;
%                     index = index + 1;
%                 end
%             end
%         end
%     end
     
    for rr = 2:numJLFsums
        r_vect = parset.powersJLF(rr,:);
        convIndices{rr}.alpha_i = [];
        %index = 1;
        for ii = 1:numAlphas
            i_vect = parset.powersMeasFunc(ii,:); %vector of exponents in the monomial
            for jj = ii:numAlphas
                j_vect = parset.powersMeasFunc(jj,:); %vector of exponents in the monomial
                if ((i_vect+j_vect)==r_vect)
                    if (ii == jj)
                        convIndices{rr}.alpha_i = [convIndices{rr}.alpha_i ii];
                    else
                        convIndices{rr}.alpha_i = [convIndices{rr}.alpha_i ii jj];
                    end
        %           index = index+1;
                end
            end
        end
        
        convIndices{rr}.alpha_i = sort(convIndices{rr}.alpha_i);
        convIndices{rr}.alpha_j = fliplr(convIndices{rr}.alpha_i);
    end
    
end

%Calculate the actual "convolution":
beta = zeros(1,numJLFsums-1);

for rr = 2:numJLFsums
    beta(rr-1) = sum(alpha(convIndices{rr}.alpha_i) .* alpha(convIndices{rr}.alpha_j));
end

end
