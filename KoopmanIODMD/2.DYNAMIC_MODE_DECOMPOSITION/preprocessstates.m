function [states,meansteadystate,scalingfactor] = preprocessstates(states)
% preprocessstates(states)

%find mean of steady state assuming stady state is from 1 to 5 time step
steadystates  = states(:,1:5);
meansteadystate = mean(steadystates,2);  

%subtract the steady state dynamics from states
states = states - meansteadystate;

% states=detrend(states);
% scalingfactor=var(var(states));
scalingfactor = 1;
states=states./scalingfactor;

states=states.^2;
end

