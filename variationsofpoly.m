%% Fall 2020
iterations = 1000000;
%%%%%%%%%%%%Starting Data%%%%%%%%%%%%%%%%
%%%Data about population
initialinds = 500; %Our population will have space for N*2 individuals
genes = 3;
traits = 3;
%%%Space for individuals 
contributions = NaN(traits,genes,initialinds*3);
alleles = NaN(initialinds*3,genes);
%%%Starting everyone the same:
%beginalleles = [1 0 0];
begincontri = [1 0 0; 0 1 0; 0 0 1]; %2
%begincontri = [1 0 0;1 0 0;1 0 0]; %1
%begincontri = [0.33 0.33 0.33;0.33 0.33 0.33; 0.33 0.33 0.33]; %3
clones = 6;
%%Gaussian Parameters:
%height = 1;
%parama = 1/2;
%paramb = 1/2; 
%ellipcenter = [0, 0]; 

%%Quadratic or Linear Parameters:
angularcoef = [1 1 1];
degrees = 2;

nasce = 0;
morre = 0;
%fid = fopen( 'results.txt', 'wt' );
%Fi = zeros(1,iterations);
%Pol = zeros(1,iterations);
%pheno = zeros(2,iterations);
%%Simulation
Fi = 0; 
Pol = 0;
Mon = 0; 
for k = 1:clones
      
     for j = 1:initialinds   
          %matrix which line is an individual    
       % alleles(j,:) = beginalleles;
        alleles(j,:) = rand(1,genes); %everyone starts with different allele matrix
        %tensor which depth is individual
        contributions(:,:,j) = begincontri;
        %contributions(:,:,k) = rand(traits,genes); 
        %rowsum = sum(contributions(:,:,k),2);
        %contributions(:,:,k) = bsxfun(@rdivide, contributions(:,:,k), rowsum);
     end
     %%Starting data
    phenotypes = phecalculator(contributions,alleles);
    fitness = fitcalculator(phenotypes, angularcoef, degrees);
    dyingvector = penalti(contributions);
    probabilities = ALLprobcalculator(contributions,alleles,angularcoef, degrees);
    [Fi,Pol] = simulate(Fi,Pol,k,iterations,fitness,contributions, probabilities, alleles, dyingvector, angularcoef, degrees); 

    %fprintf(fileID,'%f %f\n',Fi,Mon);
    

end

function [Fi,Pol] = simulate(Fi,Pol,k,iterations,fitness,contributions, probabilities, alleles, dyingvector, angularcoef, degrees) 

 for i = 1:iterations
    N = length(probabilities); 
    dB = choosing(probabilities);
    if isnan(probabilities(dB))        
            continue          
    end
    if dB <= N/2     
             [contributions,alleles,probabilities,fitness,dyingvector] = birth(dyingvector,contributions,alleles,fitness,dB, angularcoef, degrees);  
    else     
            [contributions,alleles,probabilities,fitness,dyingvector] = death(contributions,alleles,probabilities,dB,fitness,dyingvector);      
    end
    Fi(k,i) = nanmean(fitness); 
    Pol(k,i) = nanmean(sum(vecnorm(contributions,2,2))); 
   % al(i) = nanmean(vecnorm(alle));
   % alarray = [alarray, alleles];
        
    %writematrix(M, "results.txt");
    %pheno(:,i) = nanmean(vecnorm(contributions,2,2),3);  
 end

end


%%Functions



%%1. Calculating data for All Individuals
%Phenotypes
function phenotypes = phecalculator(contributions, alleles)
    %alleles are transposed
    transpose = [0,1];
    phenotypes = tmult(contributions, alleles, transpose);
    phenotypes = phenotypes(:,1,:);
    [ind, genes] = size(alleles);
    phenotypes = reshape(phenotypes,[genes,ind]);
    %~isnan(phenotypes)
    %phenotypes = phenotypes(:,1,:);
    %phenotypes = squeeze(phenotypes) %if more than 1 gene 1 trait, individuals are columns
end
  
%Fitness
function fitness = fitcalculator(phenotypes, angularcoef, degrees)  
   fitness = angularcoef*(phenotypes.^degrees) + 10;  
end


\\
%Cost - prob of dying
function dyingvector = penalti(contributions)
    %how bad is it?
    %[t,~,~] = size(contributions);
    
    strength = 0;
    N = sum(~isnan(contributions(1,1,:)));
    
    pen = 0.01 * (1 + strength*(sum(vecnorm(contributions,2,2))));
    %pen = 0.01 * (1 + strength*(t - sum(vecnorm(contributions,2,2))));
    dyingvector = reshape(pen,[],size(pen,2),1)';
    dyingvector = dyingvector*N;
end




%Creating All the probabilities 
function [probabilities,fitness,dyingvector] = ALLprobcalculator(contributions,alleles, angularcoef, degrees)
    [nind,~] = size(alleles);
    probabilities = NaN(nind,2);
    phenotypes = phecalculator(contributions,alleles);
    %2nd collumn is death
    dyingvector = penalti(contributions);
    %1st collumn is reproduction
    fitness = fitcalculator(phenotypes, angularcoef, degrees);
    ototal = nansum(fitness) + nansum(dyingvector);
    dyingvector = (dyingvector+1)/ototal;
    probbirth = (fitness+1)/ototal;
    probabilities(:,2) = dyingvector;
    probabilities(:,1) = probbirth;
    probabilities = probabilities(:)';
    %probabilities(~isnan(probabilities));
    %average = mean(fitness);
    %Scale probabilities so that its range is in the interval [0,1].
end





%Tensor multiplication
function A = tmult(A, B, transpose)
    szB = [size(B) 1];
    szA = [size(A) 1];
    if nargin < 3
        transpose = 0;
    else
        transpose = [transpose(:); 0];
        transpose = [1 2] * (transpose(1:2) ~= 0);
        if transpose == 3
            % Permutation required. Choose the matrix which permutes fastest.
            pa = any(szA(1:2) == 1);
            pb = any(szB(1:2) == 1);
            if pa || (~pb && numel(A) < numel(B))
                if ~pa
                    p = 1:numel(szA);
                    p(1:2) = p([2 1]);
                    A = permute(A, p);
                end
                szA(1:2) = szA([2 1]);
                transpose = 2;
            else
                if ~pb
                    p = 1:numel(szB);
                    p(1:2) = p([2 1]);
                    B = permute(B, p);
                end
                szB(1:2) = szB([2 1]); 
                transpose = 1;    
            end  
        end  
    end  
    switch transpose  
        case 0
            % No transposes
            A = reshape(A, szA([1:2 end 3:end-1]));
            B = reshape(B, szB([end 1:end-1]));
            dim = 2;
            szB(1) = szA(1);
        case 1
            % First matrix transposed       
            A = reshape(A, szA([1:2 end 3:end-1]));    
            B = reshape(B, szB([1 end 2:end])); 
            dim = 1;
            szB(1) = szA(2);   
        case 2 
            % Second matrix transposed    
            A = reshape(A, szA([1 end 2:end]));
            B = reshape(B, szB([end 1:end-1]));
            dim = 3;          
            szB(2) = szB(1);          
            szB(1) = szA(1);          
    end
    % Compute the output            
    A = sum(bsxfun(@times, A, B), dim);        
    % Reshape to expected size    
    clszA = [szA ones(1, numel(szB)-numel(szA))];           
    szB = [szB ones(1, numel(szA)-numel(szB))];    
    szB(3:end) = max(szB(3:end), szA(3:end));
    A = reshape(A, szB); 
end







%2. Single individual functions
%Fitness
function onefitness = onefit(singlephenotype, angularcoef, degrees) 
    onefitness = angularcoef*(singlephenotype.^degrees) + 10;
end





%%3.Event-type functions

%Choosing an Event
function chosen = choosing(probabilities)
    choices = probabilities(~isnan(probabilities)); 
    n = length(choices); 
    if any(choices)      
        dB = randsample(n,1,true,choices);     
    else  
        dB = randsample(n,1,true);     
    end
    Truenan=(isnan(probabilities)); 
    sumN=cumsum(Truenan);
    numberofnan=sumN(~isnan(probabilities));
    chosen = numberofnan(dB) + dB;
    %la= length(probabilities);
    %nanmean(probabilities(1:la/2)) >= probabilities(chosen)
end






%Death Event
 function [contributions,alleles,probabilities,fitness,dyingvector] = death(contributions,alleles,probabilities,dB,fitness,dyingvector)
    [side1,side2,~] = size(contributions);
    [~,N] = size(probabilities);
    probabilities(dB) = NaN;
    dB = dB - N/2;
    contributions(:,:,dB) = NaN(side1,side2);
    alleles(dB,:) = NaN(side2,1);
    probabilities(dB) = NaN;
    fitness(dB) = NaN;
    dyingvector(dB) = NaN;
    %recalculaitng death prob
    dyingvector = penalti(contributions);
    ototal = nansum(fitness) + nansum(dyingvector);
    dyingvector = (dyingvector +1)/ototal;
    probbirth = (fitness +1)/ototal;
    probs(:,2) = dyingvector;
    probs(:,1) = probbirth;
    probabilities = probs(:)';
 end




%Birth Event
function [contributions,alleles,probabilities,fitness,dyingvector] = birth(dyingvector, contributions,alleles, fitness,dB, angularcoef, degrees)
    cont = contributions(:,:,dB);
    alle = alleles(dB,:);
    [cont,alle] = mutate(cont,alle);
    emptyspace = find(isnan(squeeze(contributions(1,1,:))));
    emptyspace = emptyspace(1);
    contributions(:,:,emptyspace) = cont;
    alleles(emptyspace,:) = alle;
    newpheno = cont*alle';
    fitness(emptyspace) = onefit(newpheno, angularcoef, degrees);
    dyingvector = penalti(contributions);
    ototal = nansum(fitness) + nansum(dyingvector);
    dyingvector = (dyingvector+1)/ototal;
    probbirth = (fitness+1)/ototal;
    probs(:,2) = dyingvector;
    probs(:,1) = probbirth;
    probabilities = probs(:)';
end




function bool = randprob(columns, prob)
    bool = rand(columns,1) < prob;
end




function [newcont,alle] = mutate(cont,alle)
    [traits,genes,~] = size(cont);
    rate_allele_mutation = 0.05;
    badodds = 0.7;
    %mutation in the alleles
    mutation = rand(1,genes);
    whomutates = mutation<= rate_allele_mutation;
    min = -0.015;
    max = 0.015;
    if rand(1,1) > badodds
        scale_allele_mutation = (max-0).*rand(1,genes) + 0; 
    else
        scale_allele_mutation = (0-min).*rand(1,genes) + max; 
    end
    
   
    
   rate_contribution_mutation = 0.0;
    mutation = rand(traits,genes-1);
    %notmutates = ~(mutation<= rate_contribution_mutation);
    %rate_wholetrait_mutation = 0;
    %Mutation in contribution
    whomutates = mutation<= rate_contribution_mutation;
    min = -0.015;
    max = 0.015;
    scale_contribution_mutation = (max-min).*rand(traits,genes-1) + min;
    howmutates = scale_contribution_mutation.*whomutates;
    
    %cy = any(howmutates);
    %cy = any(cy);
    %cant be negative

    contcoords = cont(1:traits, 2:genes);
    
    newcontc = contcoords + howmutates;
    newcontc(newcontc < 0 ) = contcoords( newcontc < 0 );

    while sum(newcontc,2) > 1
        scale_contribution_mutation = (max-min).*rand(traits,genes-1) + min;
        howmutates = scale_contribution_mutation.*whomutates;
        newcontc = contcoords + howmutates;
        newcontc(newcontc < 0 ) = contcoords( newcontc < 0 );
    end

    newcont = [1 - sum(newcontc,2), newcontc ]; 
    
end



function X = randsphere(m,n,r)

% This function returns an m by n array, X, in which 
% each of the m rows has the n Cartesian coordinates 
% of a random point uniformly-distributed over the 
% interior of an n-dimensional hypersphere with 
% radius r and center at the origin.  The function 
% 'randn' is initially used to generate m sets of n 
% random variables with independent multivariate 
% normal distribution, with mean 0 and variance 1.
% Then the incomplete gamma function, 'gammainc', 
% is used to map these points radially to fit in the 
% hypersphere of finite radius r with a uniform % spatial distribution.
% Roger Stafford - 12/23/05

X = randn(m,n);
s2 = sum(X.^2,2);
X = X.*repmat(r*(gammainc(s2/2,n/2).^(1/n))./sqrt(s2),1,n);

end
