%SAME AS SIGMOID.m but forced ASYM = 1
%b1 old asympt removed
%b1 is E50 (or inflection pt or threshold)--fepsp at which pspike is 1/2max 
%b2 is slope (span or horizontal stretch)

function S = SIGMOID_ASYM1(b,x)
	b1 = b(1);
    b2 = b(2);
	S = 1 ./ (1 + exp((b1 - x)/b2)); 