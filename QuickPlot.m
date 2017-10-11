function s = QuickPlot(file); 

setenv("GNUTERM","qt");
a = load(file);
t = a(:,1);
w = a(:,2);

[freq,amp] = getFFT(t,w);

subplot(1,2,1);
plot(t,w);
xlim([min(t) max(t)]);
subplot(1,2,2);
loglog(freq,amp);
while (true)
pause();
end
EOF
