# Create and work with timestamps

Create and work with timestamps

## Usage

``` r
FormatTime(time, locale = TRUE, tz = "UTC", format = TSFormats(type = "R"))

Timestamp(tz = "UTC", format = TSFormats(type = "R"))

TSFormats(type = c("R", "loom"))
```

## Arguments

- time:

  A character timestamp

- locale:

  Change the timestamp of to the timezone of the locale as determined by
  [`Sys.timezone`](https://rdrr.io/r/base/timezones.html)

- tz:

  A character string specifying the time zone to be used for the
  conversion. System-specific (see
  [`as.POSIXlt`](https://rdrr.io/r/base/as.POSIXlt.html)), but `""` is
  the current time zone, and `"GMT"` is UTC. Invalid values are most
  commonly treated as UTC, on some platforms with a warning.

- format:

  A character string. The default for the `format` methods is
  `"%Y-%m-%d %H:%M:%S"` if any element has a time component which is not
  midnight, and `"%Y-%m-%d"` otherwise. If
  [`options`](https://rdrr.io/r/base/options.html)`("digits.secs")` is
  set, up to the specified number of digits will be printed for seconds.

- type:

  Type of format to get, currently supports the following formats:

  - “R”

  - “loom”

  See **Timestamp formats** below for more details

## Value

`FormatTime`: A character with the formated timestamp

`Timestamp`:A character with the current time in the specified format

`TSFormats`: Format for a specified type

## Timestamp formats

The following formats can be provided by `TSFormats`:

### R

A 24-hour R-friendly format; stores date/time information as the
following:

- Four-digit year (eg. “2020” for the year 2020)

- Two-digit month (eg. “04” for the month of April)

- Two-digit date (eg. “02” for the second day of the month)

- The letter “T”

- 24-hour two-digit time (eg. “15” for 3:00 PM)

- Two-digit minute (eg “05” for five minutes past the hour)

- Two-digit second (eg. “05” for five seconds into the minute)

- The letter “Z”

This results in a timestamp format of “YYYYMMDDTHHMMSS.SSSSSSZ” (eg.
“20200402T150505Z” for April 4th, 2020 at 3:05:05 PM). **Note**: this is
considered “R-friendly” as it does *not* use precise values for seconds

### loom

The standard format for loom timestamps; stores date/time information as
the following:

- Four-digit year (eg. “2020” for the year 2020)

- Two-digit month (eg. “04” for the month of April)

- Two-digit date (eg. “02” for the second day of the month)

- The letter “T”

- 24-hour two-digit time (eg. “15” for 3:00 PM)

- Two-digit minute (eg “05” for five minutes past the hour)

- Seconds precise to the millionth (six digits after the decimal)

- The letter “Z”

This results in a timestamp format of “YYYYMMDDTHHMMSS.SSSSSSZ” (eg.
“20200402T150505Z” for April 4th, 2020 at 3:05:05 PM). **Note**: this is
*not* considered “R-friendly” as it contains precise values for seconds.
To properly format this time for pretty-printing, please remember to
strip the precise seconds

## Examples

``` r
# \donttest{
# Get a timestamp format
srtdisk:::TSFormats()
#> [1] "%Y%m%dT%H%M%SZ"

# Create a timestamp
srtdisk:::Timestamp()
#> [1] "20260122T181537Z"

# Format a timestamp for easy viewing
time <- "20200804T214148Z"
srtdisk:::FormatTime(time)
#> [1] "2020-08-04 17:41:48 EDT"
# }
```
