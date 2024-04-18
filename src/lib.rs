use pyo3::prelude::*;
use pyo3::types::PyString;

/// Formats the sum of two numbers as string.
#[pyfunction]
#[pyo3(name = "format_read")]
fn format_read(
    _py: Python,
    title: &Bound<'_, PyString>,
    sequence: &Bound<'_, PyString>,
    quality: &Bound<'_, PyString>,
) -> PyResult<String> {
    Ok(format!("@{title}\n{sequence}\n+\n{quality}\n").to_string())
}

// #[pyfunction]
// #[pyo3(name = "asap_to_kite")]
// fn asap_to_kite(_py: Python) ->PyResult<String> {
//     Ok(())
// }

#[pyfunction]
#[pyo3(name = "rearrange_reads")]
fn rearrange_reads(
    _py: Python,
    sequence1: &Bound<'_, PyString>,
    sequence2: &Bound<'_, PyString>,
    sequence3: &Bound<'_, PyString>,
    quality1: &Bound<'_, PyString>,
    quality2: &Bound<'_, PyString>,
    quality3: &Bound<'_, PyString>,
    conjugation: &Bound<'_, PyString>,
) -> PyResult<Vec<String>> {
    let sequence1 = sequence1.to_string();
    let sequence2 = sequence2.to_string();
    let sequence3 = sequence3.to_string();

    let quality1 = quality1.to_string();
    let quality2 = quality2.to_string();
    let quality3 = quality3.to_string();

    let new_sequence1: String;
    let new_sequence2: String;
    let new_quality1: String;
    let new_quality2: String;
    let conjugation: &str = conjugation.to_str()?;

    match conjugation {
        "TotalSeqA" => {
            new_sequence1 = sequence2 + &sequence1[..10].to_string();
            new_sequence2 = sequence3;

            new_quality1 = quality2 + &quality1[..10].to_string();
            new_quality2 = quality3;
        }
        "TotalSeqB" => {
            new_sequence1 =
                sequence2 + &sequence3[..10].to_string() + &sequence3[25..34].to_string();
            new_sequence2 = sequence3[10..25].to_string();

            new_quality1 = quality2 + &quality3[..10].to_string() + &quality3[25..34].to_string();
            new_quality2 = quality3[10..25].to_string();
        }
        &_ => todo!(),
    }

    let return_values: Vec<String> = vec![new_sequence1, new_sequence2, new_quality1, new_quality2];
    Ok(return_values)
}

/// A Python module implemented in Rust.
#[pymodule]
fn asap_o_matic(_py: Python, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(format_read, m)?)?;
    m.add_function(wrap_pyfunction!(rearrange_reads, m)?)?;
    // m.add_function(wrap_pyfunction!(asap_to_kite, m)?)?;
    Ok(())
}
