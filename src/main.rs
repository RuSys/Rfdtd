extern crate Rfdtd;

use Rfdtd::fdtd2d::*;
use std::time::Instant;
use std::fs::OpenOptions;

use std::io::*;

// 待ち
pub fn waitkey() {
    let mut __ret=String::new();
    std::io::stdin().read_line(&mut __ret).ok();
    return;
}
    

#[allow(unused_must_use)]
fn main() {

    let mut file = OpenOptions::new()
            .write(true)
            .append(true)
            .create(true)
            .open("foo.txt").unwrap();

    //let mut fmodule = fdtdBuilder::newSize(1000,1000);    // セル数を指定し設定して生成
    let mut fmodule = fdtdBuilder::new();                   // デフォルトの大きさで設定して生成

    let mut t = fmodule.dt;

    fmodule.epsmu(68-20,68+20, 68-20,68+20,3.0);

    fmodule.init_source(68-40,68);

    fmodule.setup();

    let start = Instant::now();

    for _s in 1..NSTEP as usize{
        println!("Time step:{} --- Time:{}ms", _s, t);
        fmodule.e_cal();
        fmodule.feed(68-40,68, t.clone());
        fmodule.e_pml();

        t = t + 0.5 * fmodule.dt;

        fmodule.h_cal();
        fmodule.h_pml();

        t = t + 0.5 * fmodule.dt;
        
        if _s % 500 == 0 {
            write!(file, "Time step:{} --- Time:{}ms\r\n", _s, t);
            fmodule.out_file(&mut file,68-40, 68);
        }
    }

    let end = start.elapsed();
    println!("結果: {}.{:03}秒",end.as_secs(), end.subsec_nanos() / 1_000_000);

    waitkey()
}
