APID_level3_validationGOBP <- read_csv("~/Documents/Thesis/APID_level3_validationGOBP.csv")
ggplotly(ggplot(data = APID_level3_validationGOBP,
       aes(x = Iteration, y = pvalue, colour = A)) +
  geom_line() +
  scale_y_continuous(trans = "reverse") +
  labs(x = "Iteration", y = "Adjusted NEAT p-value") +
  scale_colour_discrete(name = "GO Biological\nProcess term") +
  geom_hline(aes(yintercept = 0.05),
             color = "tomato",
             linetype = "dashed",
             size = 1,
             show.legend = FALSE) +
  theme(panel.grid.minor = element_blank()))


RegNet_validationGOBP <- read_csv("~/Documents/Thesis/RegNet_validationGOBP.csv")
ggplotly(ggplot(data = RegNet_validationGOBP,
                aes(x = Iteration, y = pvalue, colour = A)) +
           geom_line() +
           scale_y_continuous(trans = "reverse") +
           labs(x = "Iteration", y = "Adjusted NEAT p-value") +
           scale_colour_discrete(name = "GO Biological\nProcess term") +
           geom_hline(aes(yintercept = 0.05),
                      color = "tomato",
                      linetype = "dashed",
                      size = 1,
                      show.legend = FALSE) +
           theme(panel.grid.minor = element_blank()))

mwtabosignal_validationGOBP <- read_csv("~/Documents/Thesis/mwtabosignal_validationGOBP.csv")
ggplotly(ggplot(data = mwtabosignal_validationGOBP,
                aes(x = Iteration, y = pvalue, colour = A)) +
           geom_line() +
           scale_y_continuous(trans = "reverse") +
           labs(x = "Iteration", y = "Adjusted NEAT p-value") +
           scale_colour_discrete(name = "GO Biological\nProcess term") +
           geom_hline(aes(yintercept = 0.05),
                      color = "tomato",
                      linetype = "dashed",
                      size = 1,
                      show.legend = FALSE) +
           theme(panel.grid.minor = element_blank()))

brain_adult_greater_mean_validationGOBP <- read_csv("~/Documents/Thesis/brain_adult_greater_mean_validationGOBP.csv")
ggplotly(ggplot(data = brain_adult_greater_mean_validationGOBP,
                aes(x = Iteration, y = pvalue, colour = A)) +
           geom_line() +
           scale_y_continuous(trans = "reverse") +
           labs(x = "Iteration", y = "Adjusted NEAT p-value") +
           scale_colour_discrete(name = "GO Biological\nProcess term") +
           geom_hline(aes(yintercept = 0.05),
                      color = "tomato",
                      linetype = "dashed",
                      size = 1,
                      show.legend = FALSE) +
           theme(panel.grid.minor = element_blank()))


brain_adult_2_5_validationGOBP <- read_csv("~/Documents/Thesis/brain_adult_2.5_validationGOBP.csv")
ggplotly(ggplot(data = brain_adult_2_5_validationGOBP,
                aes(x = Iteration, y = pvalue, colour = A)) +
           geom_line() +
           scale_y_continuous(trans = "reverse") +
           labs(x = "Iteration", y = "Adjusted NEAT p-value") +
           scale_colour_discrete(name = "GO Biological\nProcess term") +
           geom_hline(aes(yintercept = 0.05),
                      color = "tomato",
                      linetype = "dashed",
                      size = 1,
                      show.legend = FALSE) +
           theme(panel.grid.minor = element_blank()))

