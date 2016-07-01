
%save('zebtemp', 'X', 'Y') (x=IVs; y=memory)
[COEFS, SCORES, latent] = princomp(X); % find what the principal components are - a mixture of your original raw IVs
figure, plot(cumsum(latent)/sum(latent)) % how much variance explained by n components (no particular order) - allows you to see how many components you are looking for


% looking at COEFS tells you how much each of your original raw IVs loads
% onto each of the components. from the relative loading etc, you then need
% to decide what the psychological character of the component is.
% Components are order according to increasing amount of variance accounted
% for individually (1st component = most variance)

data = SCORES(:,1:5); % gives you IVs, of each of the different principal components. based on the cumsum (above), decide how many components to include - in this case 5. 
figure, 
subplot(1,2,1), scatter(SCORES(:,1), SCORES(:,2), 'or')
subplot(1,2,2), scatter3(SCORES(:,3), SCORES(:,4), SCORES(:,5), 'ob')

[b,dev,stats] = glmfit(data, Y);

figure, bar(COEFS(:,4))
