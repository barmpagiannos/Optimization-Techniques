
% Τεχνικές Βελτιστοποίησης
% 2η εργαστηριακή άσκηση
% Μπαρμπαγιάννος Βασίλειος
% ΑΕΜ: 10685

clc, clear, close all;

% Ορισμός πεδίου ορισμού.
% Δημιουργούμε ένα πλέγμα (grid) σημείων (x,y) για να σχεδιάσουμε τη συνάρτηση.
% Επιλέγουμε το διάστημα [-3, 3] που καλύπτει τα σημεία ενδιαφέροντος.
% (0,0), (1,1), (-1,-1) που αναφέρει η εκφώνηση.
x = linspace(-2.5, 2.5, 100); % 100 σημεία από το -2.5 έως το 2.5.
y = linspace(-2.5, 2.5, 100);
[X, Y] = meshgrid(x, y); % Δημιουργία του 2D πλέγματος πάνω στο οποίο θα σχεδιαστεί η συνάρτηση.

% Υπολογισμός της συνάρτησης f(x,y).
% f(x,y) = x^3 * exp(-x^2 - y^4)
Z = (X.^3) .* exp(-X.^2 - Y.^4);

% Γραφική παράσταση 3D (Surface Plot).
figure('Name', 'Theme 1: 3D Surface Plot', 'NumberTitle', 'off');
surf(X, Y, Z); 
shading interp;
colormap jet;
colorbar;
xlabel('x');
ylabel('y');
zlabel('f(x,y)');
title('3D Γράφημα της f(x,y) = x^3 e^{-x^2 - y^4}');
grid on;

% Γραφική παράσταση ισοϋψών καμπυλών (Contour Plot).
% Αυτό το γράφημα δείχνει τις ισοϋψείς καμπύλες.
figure('Name', 'Theme 1: Contour Plot', 'NumberTitle', 'off');
contourf(X, Y, Z, 20, 'HandleVisibility', 'off'); % 20 ισοϋψείς καμπύλες.
colorbar;
xlabel('x');
ylabel('y');
title('Ισοϋψείς Καμπύλες της f(x,y)');
grid on;

% Προσθήκη των σημείων εκκίνησης που θα χρησιμοποιήσουμε αργότερα.
hold on;
plot(0, 0, 'ro', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '(0,0)');
plot(-1, -1, 'gx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '(-1,-1)');
plot(1, 1, 'mx', 'MarkerSize', 10, 'LineWidth', 2, 'DisplayName', '(1,1)');
legend show;
hold off;
