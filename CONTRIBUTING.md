# Contributing

Contributions are welcome, please open a pull request with any changes and add Dialpuri as a reviewer. 

## Development

To build NucleoFind upon every Python import, which will be useful in development, run:

    pip install --no-build-isolation --config-settings=editable.rebuild=true -Cbuild-dir=build -ve .


    
## Testing
Any changes must pass the tests defined in `package/tests`. Test can be ran using `pytest` with: 

    pytest package/tests --unit --runslow -v

