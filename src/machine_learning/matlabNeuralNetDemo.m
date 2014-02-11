% Basic (Matlab) neural net demo

% Graph logsig function:
n = -5:0.1:5;
a = logsig(n);
%figure; plot(n,a)

% Input features, just integers
P = [0 1 2 3 4 5 6 7 8 9 10];

% T = target output
T = [0 1 2 3 4 3 2 1 2 3 4];

% Create a new network: 
%    range of input is [0,10]
%    has 5 neurons in a single hidden layer
%    and a single output neuron
%    transfer functions of hidden layer are logsig: 
%        logsig(n) = 1 / (1 + exp(-n))
% Could also use tansig: tansig(n) = 2/(1+exp(-2*n))-1 
%    transfer functions of output layer are purelin
net = newff([0 10],[5 1],{'logsig' 'purelin'});

% Show the output before training.
% sim runs the network called net on data set P.
Y = sim(net,P);
figure; hold on;
hTarget = plot(P,T);
hNet = plot(P,Y,'ro');
title('Untrained Network');
legend([hTarget, hNet], 'Target', 'Output of network]');
hold off;

% Train the network for a specified time
net.trainParam.epochs = 50;
net = train(net,P,T);

% Run the network on the input again to see how well it learned it.
Y = sim(net,P);
figure; hold on;
hTarget = plot(P,T);
hNet = plot(P,Y,'ro');
hold off;


% What's really going on for other data points? 
% What function did it learn?
P2 = [0:0.1:10];
Y2 = sim(net,P2);
figure; hold on;
hTarget = plot(P,T);
hNet = plot(P2,Y2,'ko');
hold off;

title('Trained Network');
legend([hTarget, hNet], 'Target', 'Output of network]');
hold off;

% net.IW{1} ({} since IW is a cell array)
% net.LW{2,1}

% Input weights:
%    -1.0156
%    -1.3557
%    -0.8843
%    -1.0095
%    -7.8035
% 
% Input biases:
%     8.0644
%   -21.9935
%     4.8510
%     3.1040
%    -1.2886
% 
% Layer weights:
%     -8.0044    
%     1.3550   
%     11.8593   
%     -8.1436   
%     -2.8339
%    
% Output bias:
%     4.6408

   
   