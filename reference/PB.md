# Create a progress bar

Progress bars are useful ways of getting updates on how close a task is
to completion. However, they can get in the way of RMarkdown documents
with lots of unnecesssary printing. `PB` is a convenience function that
creates progress bars with the following defaults

- `char = '='`

- `style = 3`

- `file = stderr()`

## Usage

``` r
PB()
```

## Value

An object of class
[`txtProgressBar`](https://rdrr.io/r/utils/txtProgressBar.html)

## See also

[`txtProgressBar`](https://rdrr.io/r/utils/txtProgressBar.html)
[`stderr`](https://rdrr.io/r/base/showConnections.html)

## Examples

``` r
# \donttest{
pb <- srtdisk:::PB()
for (i in 1:10) {
  utils::setTxtProgressBar(pb, i / 10)
}
close(pb)
# }
```
