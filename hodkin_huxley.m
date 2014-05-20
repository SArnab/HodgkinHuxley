function[timeVector,V,extI,ionNa,ionK,netIon,spikeCount] = hodkin_huxley(stimulusCurrent)
    tOnset = 40;
    tMax = 200;
    tStep = 0.01;
    
    e_K = -85;
    v_L = -65;
    e_Na = 99;
    
    g_Na = 400;
    g_K = 200;
    g_L = 2;
    c = 2;
    i0 = stimulusCurrent;
    threshold = -30;
    spikeCount = 0;
    
    timeVector = (0:tStep:tMax);
    
    V = zeros(1,length(timeVector));
    
    alpha_m = (0.1*(V(1)+40))/(1-exp(-0.1*(V(1)+40)));
    alpha_h = .07*exp(-.05*(V(1)+65));    
    alpha_n = (.01*(V(1)+55))/(1-exp(-0.1*(V(1)+55)));
    
    beta_m = 4*exp(-.0556*(V(1)+65));
    beta_h = 1/(1+exp(-0.1*(V(1)+35)));
    beta_n = .125*exp(-.0125*(V(1)+65));
    
    m = V;
    m(1) = alpha_m/(alpha_m + beta_m);
    
    h = V;
    h(1) = alpha_h/(alpha_h+beta_h);
    
    n = V;
    n(1) = alpha_n/(alpha_n + beta_n);
    
    extI = V;
    extI(timeVector >= tOnset) = i0;
    
    ionNa = V;
    ionK = V;
    ionL = V;
    netIon = V;
    
    if(V(1) > threshold)
        aboveThreshold = true;
    else
        aboveThreshold = false;
    end
    
    tic;
    for i=1:length(timeVector)-1
        if(V(i) >= threshold && aboveThreshold == false )
            %fprintf('Threshold crossed at time %d\n',timeVector(i));
            aboveThreshold = true;
            spikeCount = spikeCount + 1;
        elseif(V(i) < threshold && aboveThreshold == true)
            aboveThreshold = false;
        end
        % Calculate gating coefficients at this instant
        alpha_m = (0.1*(V(i)+40))/(1-exp(-0.1*(V(i)+40)));
        alpha_n = (.01*(V(i)+55))/(1-exp(-0.1*(V(i)+55)));
        alpha_h = .07*exp(-.05*(V(i)+65));
        beta_m = 4*exp(-.0556*(V(i)+65));
        beta_h = 1/(1+exp(-0.1*(V(i)+35)));
        beta_n = .125*exp(-.0125*(V(i)+65));
        
        % Estimate gating probabilities using Forward Euler
        m(i+1) = m(i) + tStep*((alpha_m*(1-m(i)))-(beta_m*m(i)));
        h(i+1) = h(i) + tStep*((alpha_h*(1-h(i)))-(beta_h*h(i)));
        n(i+1) = n(i) + tStep*((alpha_n*(1-n(i)))-(beta_n*n(i)));
        
        % Calculate Ion Contributions
        ionK(i) = g_K * (n(i)^4) * (V(i)-e_K);
        ionNa(i) = g_Na * (m(i)^3) * h(i) * (V(i)-e_Na);
        ionL(i) = g_L * (V(i)-v_L);
        netIon(i) = extI(i) - ionL(i) - ionNa(i) - ionK(i);
        
        % Estimate membrane potential using Forward Euler
        V(i+1) = V(i) + tStep*(netIon(i)/c);
    end
    toc;
    %plot(timeVector,V);
   % disp(spikeCount);
end