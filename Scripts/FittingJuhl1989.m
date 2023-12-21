data = readtable('../data/Juhl_ClinSci1989_fig1.csv')

plot(data.Day, data.P_PV, '*-')
xd = data.Day(2:end)
yd = data.P_PV(2:end)

     f(x) = a + b*exp(-c*x)
Coefficients (with 95% confidence bounds):
       a =       10.97  (7.349, 14.6)
       b =       9.518  (4.933, 14.1)
       c =     0.05833  (-0.01525, 0.1319)
       
       y = 10.95831 + 8.991122*e^(-0.05793506*x)