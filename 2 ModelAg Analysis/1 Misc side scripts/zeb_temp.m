
%save('zebtemp', 'X', 'Y')

[COEFS, SCORES, latent] = princomp(X)
figure, plot(cumsum(latent)/sum(latent))

data = SCORES(:,5);
figure, 
subplot(1,2,1), scatter(SCORES(:,1), SCORES(:,2), 'or')
subplot(1,2,2), scatter3(SCORES(:,3), SCORES(:,4), SCORES(:,5), 'ob')



