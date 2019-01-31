eff = seq(0.1, 1, 0.01)
r0_vals = c(1.1,1.3,1.6,2.0,2.5,4.5)
vc = NULL

hit_thresh = function(r0, eff){
  vc = (1 - (1/r0))/eff
}

saved = matrix(0, nrow = length(r0_vals), ncol = length(eff))

for (i in 1:length(r0_vals)){
  for (j in 1:length(eff)){
    saved[i, j] = (1 - (1/r0_vals[i]))/eff[j]
  }
}

plot(eff, saved[1,], type = "l", ylim = c(0,1), xlim = c(0.1,1.1), xaxt = "n", yaxt = "n",
     xlab = "Vaccine effectiveness", ylab = "Critical vaccination coverage")
axis(side = 1, at = seq(0.1,1,0.1))
axis(side = 2, at = seq(0, 1, 0.1))
text(x = 1, y = saved[1,91], labels = bquote("R"[0] == 1.1), pos = 4)
lines(eff, saved[2,])
text(x = 1, y = saved[2,91], labels = bquote("R"[0] == 1.3), pos = 4)
lines(eff, saved[3,])
text(x = 1, y = saved[3,91], labels = bquote("R"[0] == 1.6), pos = 4)
lines(eff, saved[4,])
text(x = 1, y = saved[4,91], labels = bquote("R"[0] == 2.0), pos = 4)
lines(eff, saved[5,])
text(x = 1, y = saved[5,91], labels = bquote("R"[0] == 2.5), pos = 4)
lines(eff, saved[6,])
text(x = 1, y = saved[6,91], labels = bquote("R"[0] == 4.5), pos = 4)


library(reshape2)

df = data.frame(r0_vals, saved)
df_melted = melt(df, id.vars = "r0_vals")

ggplot(df_melted, aes(x = variable, y = value)) + 
  geom_line(aes(color = r0_vals, group = r0_vals)) + 
  scale_x_continuous(breaks = eff)

df = data.frame(cat = LETTERS[1:6], VAR1 = runif(6), VAR2 = runif(6), VAR3 = runif(6), VAR4 = runif(6))
df_melted = melt(df, id.vars = 'cat')

