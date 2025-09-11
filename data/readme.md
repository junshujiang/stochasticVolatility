



## Daily Data

### IF index
Range 2017-2024, data source is from JBES paper.

### Option for IF index

Option chain price observed at 20240417. `calls_records.csv` and `put_records.csv`


`IFtest.csv` is the data between `2024-03-25` are `2024-04-19`.

price of SPY option and its underlying at 9th Sep, 2025.


### Option
```MySql
Run query SELECT * FROM `option_chain` WHERE `act_symbol` = 'SPY'
```

preprocessed data: `SPY_testPeriod.csv` selects one months


### SPY underlying assets
```MySql
Run query SELECT * FROM `ohlcv` WHERE `act_symbol`='SPY'
```
`SPY.csv`

`SPV_option.csv`


## HFT data


### IF HFT 
`weighted_priceIF2404_0317_417` 0315,2024 - 0415,2024