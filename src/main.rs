use ark_ff::MontFp;
use ark_poly::univariate::DensePolynomial;
//pub mod field;
pub mod sidh_sike_p434;
use sidh_sike_p434::{Fq2 as F, F as Fp};
pub mod matrix;
use matrix::*;
use merkle::{poseidon_parameters, FieldMT, FieldPath};
use std::{
    fs::{self, File},
    io::{self, Write},
    path::Path,
    str::FromStr,
    time::Instant,
};

// TODO: Move to separate crate
pub mod generalized_fri;
pub mod get_roots;
pub mod isogeny_prove;
use ark_crypto_primitives::crh::poseidon;
use ark_crypto_primitives::CRHScheme;
use ark_ff::Field;
use ark_ff::UniformRand;
use ark_poly::Polynomial;
use ark_std::test_rng;
use generalized_fri::*;
use isogeny_prove::{prove, verify};
pub mod merkle;
use ark_serialize::{CanonicalSerialize, Compress};
fn main() -> io::Result<()> {
    let mut results = Vec::new();
    for i in 0..5 {
        if let Some(result) = round(i) {
            results.push((i, result.0, result.1, result.2));
        }
    } // Specify the path to the file where you want to save the results
    let file_path = "results_new.txt";
    write_results_to_file(&results, file_path)?;

    Ok(())
}

fn lines_from_file(filename: impl AsRef<Path>) -> io::Result<Vec<F>> {
    let file_content = fs::read_to_string(filename)?;
    let lines = file_content.split(',').collect::<Vec<_>>();
    let mut results = Vec::new();
    for line in lines {
        let a: Fp;
        let b: Fp;
        if !line.contains("*a") {
            println!("line: {:?}", line);
            a = Fp::from_str(line).map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse Fp4"))?;
            b = Fp::from(0);
        } else if !line.contains("+") {
            let mut parts = line.trim().split("*a");
            b = Fp::from_str(parts.next().unwrap().trim())
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse Fp2"))?;
            a = Fp::from(0);
        } else {
            let mut parts = line.trim().split("*a +");
            b = Fp::from_str(parts.next().unwrap())
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse Fp1"))?;
            a = Fp::from_str(parts.next().unwrap().trim())
                .map_err(|_| io::Error::new(io::ErrorKind::InvalidData, "Failed to parse Fp5"))?;
        }
        results.push(F::new(a, b));
    }
    Ok(results)
}
fn write_results_to_file(results: &Vec<(usize, u64, u64, f32)>, file_path: &str) -> io::Result<()> {
    let mut file = File::create(file_path)?;
    for result in results {
        writeln!(file, "{},{},{},{}", result.0, result.1, result.2, result.3)?;
    }
    Ok(())
}
fn round(i: usize) -> Option<(u64, u64, f32)> {
    let now = Instant::now();
    let l_list: Vec<usize> = vec![vec![2; 11 + i]].concat();
    // make s such that s_order == 2^5 times bigger than poly
    let mut s = F::new(MontFp!("20166910023067061949242007329499498203359431897882042367010355835922513064073298943410517857055490103593376746276366289375459766625"), MontFp!("9070588010778744328358501144613351095932673883221130019470787206592208103281447479286045237519441445473142536447684307414402867833"));

    for _ in 0..5 - i {
        s = s.pow(&[2]);
    }
    let mut g: F = s.clone();
    for _ in 0..5 {
        g = g.pow(&[2]);
    }
    //for _ in 0..2 {
    //    g = g.pow(&[3]);
    //}
    let r: F = F::new(Fp::from(5), Fp::from(3));
    let witness: DensePolynomial<F> =
        DensePolynomial { coeffs: lines_from_file(&format!("Phis/polynomial_{}.txt", 9 + i)).unwrap() };
    let n = witness.coeffs.len();

    let mut rng = test_rng();
    let a: F = F::rand(&mut rng);
    let b: F = F::rand(&mut rng);
    let c: F = F::rand(&mut rng);
    let blinding_factor: DensePolynomial<F> = DensePolynomial { coeffs: vec![a, b, c] }.naive_mul(&DensePolynomial {
        coeffs: vec![vec![-F::from(1)], vec![F::from(0); n - 1], vec![F::from(1)]].concat(),
    });
    let b_witness: DensePolynomial<F> = witness.clone() + blinding_factor;

    // psi
    let psi: DensePolynomial<F> =
        DensePolynomial { coeffs: lines_from_file(&format!("Psis/polynomial_{}.txt", 9 + i)).unwrap() };

    let y_start: F = b_witness.evaluate(&F::from(1));
    let y_end: F = b_witness.evaluate(&g.pow(&[n as u64 - 1]));

    let s_ord: u64 = n as u64 * 32;
    let rep_param: usize = 1;
    let grinding_param: u8 = 32;

    let (challenge_vals, roots_fri, roots, paths_fri, points_fri, additional_paths_and_points, ws, paths_and_points) =
        prove(witness, psi, g, s, r, s_ord, &y_start, &y_end, l_list.clone(), rep_param, grinding_param);
    println!("Prover Time: {} s", now.elapsed().as_secs());
    let prover_time = now.elapsed().as_secs();

    let now = Instant::now();
    let b = verify(
        challenge_vals.clone(),
        roots_fri.clone(),
        roots.clone(),
        paths_fri.clone(),
        points_fri.clone(),
        additional_paths_and_points.clone(),
        ws.clone(),
        g,
        s,
        r,
        &(n as u64),
        s_ord,
        &y_start,
        &y_end,
        l_list,
        rep_param,
        grinding_param,
        paths_and_points.clone(),
    );
    println!("Verifier Time: {} ms", now.elapsed().as_millis());
    let verifier_time = now.elapsed().as_millis() as u64;
    if b {
        let size1 = challenge_vals.serialized_size(Compress::Yes);
        let size2 = roots.serialized_size(Compress::Yes);
        let size3 = points_fri.serialized_size(Compress::Yes);
        let size4 = roots_fri.serialized_size(Compress::Yes);
        let size5 = paths_fri.serialized_size(Compress::Yes);
        let size6 = additional_paths_and_points[0].paths_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[0].points_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[0].paths_minus.serialized_size(Compress::Yes)
            + additional_paths_and_points[0].points_minus.serialized_size(Compress::Yes);
        let size7 = additional_paths_and_points[1].paths_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[1].points_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[1].paths_minus.serialized_size(Compress::Yes)
            + additional_paths_and_points[1].points_minus.serialized_size(Compress::Yes);
        let size8 = additional_paths_and_points[2].paths_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[2].points_plus.serialized_size(Compress::Yes)
            + additional_paths_and_points[2].paths_minus.serialized_size(Compress::Yes)
            + additional_paths_and_points[2].points_minus.serialized_size(Compress::Yes);
        let size9 = paths_and_points[0].root_f.serialized_size(Compress::Yes)
            + paths_and_points[0].root_g.serialized_size(Compress::Yes)
            + paths_and_points[0].root_u.serialized_size(Compress::Yes)
            + paths_and_points[0].merkle_paths_f.serialized_size(Compress::Yes)
            + paths_and_points[0].merkle_paths_g.serialized_size(Compress::Yes)
            + paths_and_points[0].merkle_paths_u.serialized_size(Compress::Yes)
            + paths_and_points[0].queried_points_f.serialized_size(Compress::Yes)
            + paths_and_points[0].queried_points_u.serialized_size(Compress::Yes)
            + paths_and_points[0].queried_points_g.serialized_size(Compress::Yes);
        let size10 = paths_and_points[1].root_f.serialized_size(Compress::Yes)
            + paths_and_points[1].root_g.serialized_size(Compress::Yes)
            + paths_and_points[1].root_u.serialized_size(Compress::Yes)
            + paths_and_points[1].merkle_paths_f.serialized_size(Compress::Yes)
            + paths_and_points[1].merkle_paths_g.serialized_size(Compress::Yes)
            + paths_and_points[1].merkle_paths_u.serialized_size(Compress::Yes)
            + paths_and_points[1].queried_points_f.serialized_size(Compress::Yes)
            + paths_and_points[1].queried_points_u.serialized_size(Compress::Yes)
            + paths_and_points[1].queried_points_g.serialized_size(Compress::Yes);
        let size11 = paths_and_points[2].root_f.serialized_size(Compress::Yes)
            + paths_and_points[2].root_g.serialized_size(Compress::Yes)
            + paths_and_points[2].root_u.serialized_size(Compress::Yes)
            + paths_and_points[2].merkle_paths_f.serialized_size(Compress::Yes)
            + paths_and_points[2].merkle_paths_g.serialized_size(Compress::Yes)
            + paths_and_points[2].merkle_paths_u.serialized_size(Compress::Yes)
            + paths_and_points[2].queried_points_f.serialized_size(Compress::Yes)
            + paths_and_points[2].queried_points_u.serialized_size(Compress::Yes)
            + paths_and_points[2].queried_points_g.serialized_size(Compress::Yes);

        println!(
            "Proof Size: {} kB",
            ((size1 + size2 + size3 + size4 + size5 + size6 + size7 + size8 + size9 + size10 + size11) as f32)
                / 1000f32
        );
        let proof_size = ((size1 + size2 + size3 + size4 + size5 + size6 + size7 + size8 + size9 + size10 + size11)
            as f32)
            / 1000f32;
        println!("Verification successful");
        Some((prover_time, verifier_time, proof_size))
    } else {
        println!("Verification failed");
        None
    }
}
