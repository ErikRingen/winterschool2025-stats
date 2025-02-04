# winterschool2025-stats

## Computational Environment

To run this code, you will need to [install R](https://www.r-project.org/) and the following R packages:

```         
install.packages(
c("tidyverse",
  "marginaleffects",
  "targets",
  "knitr",
  "mgcv"
)
```

The `.qmd` file is Quarto markdown, which will render the slides from the winter school session. Quarto markdown can be rendered directly from the RStudio IDE, or from the command line.

* Download Quarto from https://quarto.org/ and install it under the 'getting started' page. There is also a link to docs on how to use Quarto with either VS Code, Jupyter, RStudio, etc.

### Render a Quarto file

If you are using RStudio or VS Code then the instructions will tell you about a preview button. Clicking this will render in your selected format and open the document either within the IDE or in your browser.

Alternatively, you can use the command line to render the document. This is done with the following command:

```bash
quarto render docs/marginaleffects.qmd
```
