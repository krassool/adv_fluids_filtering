function myprint(variables)

n_vars = length(variables);

for i=1:n_vars
   pretty(variables(i))
   fprintf('\n\n')
end