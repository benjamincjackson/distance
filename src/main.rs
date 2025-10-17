use clap::parser::ValueSource;
use distance::*;

fn main() -> Result<(), DistanceError> {
    let m = get_cli_arguments();

    if let Some(ValueSource::CommandLine) = m.value_source("licenses") {
        println!("{}", licences());
        std::process::exit(0);
    }

    let setup = set_up(&m)?;
    run(setup)?;

    Ok(())
}

fn licences() -> String {
    "
Copyright 2022, Ben Jackson

distance is licensed under the GNU LIBRARY GENERAL PUBLIC LICENSE, Version 2

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

distance incorporates some of Rust-Bio, which is licensed under the MIT licence:

The MIT License (MIT)

Copyright (c) 2016 Johannes Köster, the Rust-Bio team, Google Inc.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and 
associated documentation files (the \"Software\"), to deal in the Software without restriction, 
including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, 
and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, 
subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial 
portions of the Software.

THE SOFTWARE IS PROVIDED \"AS IS\", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT 
NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND 
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES 
OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN 
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.".to_string()
}

