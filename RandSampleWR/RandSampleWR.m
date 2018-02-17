% Random smapling with weight and without replacement 
% O/P : sample--> random sample
% I/P : N--> population; n--> sample size; W--> weightage vector

function sample = RandSampleWR(N,n,W)
    %% checking
    if(length(N)<=0)
        disp('Population size must be greater than zero.')
        close
    elseif(n<=0)
        disp('Sample size must be greater than zero.')
        close
    elseif(length(N)<=n)
        disp('Population size must begreater than sample size.')
        close
    end
    if (length(N) ~= length(W))
        disp('Population size and weightage vector size must be same.')
        close
    end
    %% initialization
    maxN = length(N);
    k = 1;
    W = W/sum(W(1:maxN));
    %% Knuth-Fisher-Yates sampling
    for i=n:-1:1
        r1 = ceil(rand*maxN);
        r2 = rand;
        while(r2>W(r1))
            r1 = ceil(rand*maxN);
            r2 = rand;
        end
        sample(k) = N(r1);
        ttemp = N(maxN);
        ttemp1 = W(r1);
        ttemp2 = W(maxN);
        N(maxN) = sample(k);
        N(r1) = ttemp;
        W(maxN) = ttemp1;
        W(r1) = ttemp2;
        maxN = maxN-1;
        W = W/sum(W(1:maxN));
        k = k+1;
        clear ttemp ttemp1 ttemp2 r1 r2
    end