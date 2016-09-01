function [ts,ns] = poisson_fixed_time ( lambdas, time )

%*****************************************************************************80
%
%% POISSON_FIXED_TIME counts the Poisson events in a fied time.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    21 November 2015
%
%  Author:
%
%    John Burkardt
%    HEAVILY modified by John Barrett to allow simulating multiple Poisson
%    processes simultaneously.
%    Originally downloaded from: https://people.sc.fsu.edu/~jburkardt/m_src/poisson_simulation/poisson_simulation.html
%
%  Parameters:
%
%    Input, real LAMBDA, the average number of events per unit time.
%
%    Input, real TIME, the amount of time to observe.
%
%    Output, integer N, the number of Poisson events observed.
%
    nl = numel(lambdas);
    sl = size(lambdas);
    
    lambdas = lambdas(:);
    
    generateSteps = @() -log(rand(nl,1))./lambdas;
    dts = generateSteps();
    timeLeft = repmat(time,nl,1) - dts;
    ts = num2cell(dts);
    ts(timeLeft < 0) = repmat({[]},sum(timeLeft < 0),1);

    while any(timeLeft > 0)
        dts = generateSteps();
        ts(dts <= timeLeft) = arrayfun(@(t,dt) [t{1}; t{1}(end)+dt],ts(dts < timeLeft),dts(dts < timeLeft),'UniformOutput',false);
        timeLeft = timeLeft - dts;
    end
    
    assert(all(cellfun(@(t) all(t <= time) & issorted(t),ts)));
    
    ts = reshape(ts,sl);

  	ns = cellfun(@(t) numel(t),ts);
end