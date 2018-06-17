/*
 * This software is released under the MIT License.
 * http://opensource.org/licenses/mit-license.php
 */


use std::fs::File;
use std::io::*;

// デフォルト解析領域分割数
const NX0: i32 = 120;
const NY0: i32 = 120;

// デフォルトセルサイズ
const DX: f32 = 0.005;
const DY: f32 = 0.005;

// 計算ステップ総数
#[allow(dead_code)]
pub const NSTEP: i32 = 2000;

// pml次数, 層数, 要求精度
const LPML: i32 = 8;
const ORDER: i32 = 4;
const RMAX: f32 = -120.0;// (dB)

const NX: i32 = NX0 + 2 * LPML;
const NY: i32 = NY0 + 2 * LPML;

const NX_P: usize = NX as usize;
const NY_P: usize = NY as usize;


// 背景媒質
#[allow(non_upper_case_globals)]
const epsbk: f32 = 1.0;
#[allow(non_upper_case_globals)]
const mubk: f32 = 1.0;
#[allow(non_upper_case_globals)]
const sigebk: f32 = 0.0;
#[allow(non_upper_case_globals)]
const sigmbk: f32 = 0.0;
#[allow(non_upper_case_globals)]
const copml: f32 = -1.5280063e-4;

// 定数
const EPS0: f32 = 8.8541878e-12;
const MU0: f32 = 1.2566371e-6;
const C: f32 = 2.9979246e8;

// PML領域の位置格納構造体
#[derive(Clone,Debug)]
#[allow(non_camel_case_types)]
struct pml {
    pub x_s: i32, // x座標開始位置
    pub x_l: i32, // x座標終端位置
    pub y_s: i32, // y座標開始位置
    pub y_l: i32, // y座標終端位置
}

#[allow(non_camel_case_types)]
pub struct fdtd {
    // 時間ステップサイズ, 時間
    pub dt: f32,

    // 電界配列
    ex: Vec<Vec<f32>>,
    ey: Vec<Vec<f32>>,
    pub ez: Vec<Vec<f32>>,

    // 磁界配列
    pub hx: Vec<Vec<f32>>,
    pub hy: Vec<Vec<f32>>,
    hz: Vec<Vec<f32>>,

    // 係数配列
    aexpml: Vec<Vec<f32>>, // PML用
    aeypml: Vec<Vec<f32>>, // PML用

    aex: Vec<Vec<f32>>,
    aey: Vec<Vec<f32>>,
    aez: Vec<Vec<f32>>,

    bexpml: Vec<Vec<f32>>, // PML用
    beypml: Vec<Vec<f32>>, // PML用

    bexy: Vec<Vec<f32>>,
    beyx: Vec<Vec<f32>>,
    bezx: Vec<Vec<f32>>,
    bezy: Vec<Vec<f32>>,

    amxpml: Vec<Vec<f32>>, // PML用
    amypml: Vec<Vec<f32>>, // PML用

    amx: Vec<Vec<f32>>,
    amy: Vec<Vec<f32>>,
    amz: Vec<Vec<f32>>,

    bmxpml: Vec<Vec<f32>>, // PML用
    bmypml: Vec<Vec<f32>>, // PML用

    bmxy: Vec<Vec<f32>>,
    bmyx: Vec<Vec<f32>>,
    bmzx: Vec<Vec<f32>>,
    bmzy: Vec<Vec<f32>>,

    expml: Vec<Vec<f32>>, // PML用
    eypml: Vec<Vec<f32>>, // PML用
    ezx: Vec<Vec<f32>>, // PML用
    ezy: Vec<Vec<f32>>, // PML用

    hxpml: Vec<Vec<f32>>, // PML用
    hypml: Vec<Vec<f32>>, // PML用
    hzx: Vec<Vec<f32>>, // PML用
    hzy: Vec<Vec<f32>>, // PML用

    // 比誘電率, 導電率
    epsd: Vec<Vec<f32>>,
    sgmed: Vec<Vec<f32>>,

    // 比透磁率, 磁気伝導率
    mud: Vec<Vec<f32>>,
    sgmmd : Vec<Vec<f32>>,

    // PML領域
    pml_s: Vec<pml>,

    // 給電係数
    befed: f32,
    duration: f32,
    t0: f32,

    // セルサイズ設定
    nx: i32,
    ny: i32,
}

#[allow(non_camel_case_types)]
pub struct fdtdBuilder;

pub trait New<T> {
    #[allow(non_camel_case_types, non_snake_case)]
    fn newSize(T, T) -> fdtd;

    fn new() -> fdtd;
}


// セルサイズ設定可能な引数
#[allow(non_snake_case)]
impl New<i32> for fdtdBuilder {
    fn newSize(x:i32, y:i32) -> fdtd {
        let nxa = (x + 2 * LPML) as usize;
        let nya = (y + 2 * LPML) as usize;

        let mut n = fdtd{dt: 0.0, ex:vec![vec![0.0; nxa]; nya], ey:vec![vec![0.0; nxa]; nya], ez:vec![vec![0.0; nxa]; nya], hx:vec![vec![0.0; nxa]; nya], hy:vec![vec![0.0; nxa]; nya], hz:vec![vec![0.0; nxa]; nya],
        aexpml: vec![vec![0.0; nxa]; nya], aeypml: vec![vec![0.0; nxa]; nya], aex: vec![vec![0.0; nxa]; nya], aey: vec![vec![0.0; nxa]; nya], aez: vec![vec![0.0; nxa]; nya],
        bexpml: vec![vec![0.0; nxa]; nya], beypml: vec![vec![0.0; nxa]; nya], bexy: vec![vec![0.0; nxa]; nya], beyx: vec![vec![0.0; nxa]; nya], bezx: vec![vec![0.0; nxa]; nya], bezy: vec![vec![0.0; nxa]; nya],
        amxpml: vec![vec![0.0; nxa]; nya], amypml: vec![vec![0.0; nxa]; nya], amx: vec![vec![0.0; nxa]; nya], amy: vec![vec![0.0; nxa]; nya], amz: vec![vec![0.0; nxa]; nya],
        bmxpml: vec![vec![0.0; nxa]; nya], bmypml: vec![vec![0.0; nxa]; nya], bmxy: vec![vec![0.0; nxa]; nya], bmyx: vec![vec![0.0; nxa]; nya], bmzx: vec![vec![0.0; nxa]; nya], bmzy: vec![vec![0.0; nxa]; nya],
        expml: vec![vec![0.0; nxa]; nya], eypml: vec![vec![0.0; nxa]; nya], ezx: vec![vec![0.0; nxa]; nya], ezy: vec![vec![0.0; nxa]; nya],
        hxpml: vec![vec![0.0; nxa]; nya], hypml: vec![vec![0.0; nxa]; nya], hzx: vec![vec![0.0; nxa]; nya], hzy: vec![vec![0.0; nxa]; nya], 
        epsd: vec![vec![epsbk; nxa+1]; nya+1], sgmed: vec![vec![sigebk; nxa+1]; nya+1],
        mud: vec![vec![mubk; nxa+1]; nya+1], sgmmd : vec![vec![sigmbk; nxa+1]; nya+1],
        pml_s: Vec::new(), befed: 0.0, duration: 0.0, t0: 0.0, nx: x + 2 * LPML, ny: y + 2 * LPML};

        println!("Set cell size x:{}, y:{}", x, y);
        //時間ステップ
        let v = C / ((epsbk * mubk).sqrt());
        n.dt = 0.99999/(v * ((1.0 / (DX * DX) + 1.0 / (DY * DY)).sqrt()));

        return n;
    }

    fn new() -> fdtd {
        let mut n = fdtd{dt: 0.0, ex:vec![vec![0.0; NX_P]; NY_P], ey:vec![vec![0.0; NX_P]; NY_P], ez:vec![vec![0.0; NX_P]; NY_P], hx:vec![vec![0.0; NX_P]; NY_P], hy:vec![vec![0.0; NX_P]; NY_P], hz:vec![vec![0.0; NX_P]; NY_P],
        aexpml: vec![vec![0.0; NX_P]; NY_P], aeypml: vec![vec![0.0; NX_P]; NY_P], aex: vec![vec![0.0; NX_P]; NY_P], aey: vec![vec![0.0; NX_P]; NY_P], aez: vec![vec![0.0; NX_P]; NY_P],
        bexpml: vec![vec![0.0; NX_P]; NY_P], beypml: vec![vec![0.0; NX_P]; NY_P], bexy: vec![vec![0.0; NX_P]; NY_P], beyx: vec![vec![0.0; NX_P]; NY_P], bezx: vec![vec![0.0; NX_P]; NY_P], bezy: vec![vec![0.0; NX_P]; NY_P],
        amxpml: vec![vec![0.0; NX_P]; NY_P], amypml: vec![vec![0.0; NX_P]; NY_P], amx: vec![vec![0.0; NX_P]; NY_P], amy: vec![vec![0.0; NX_P]; NY_P], amz: vec![vec![0.0; NX_P]; NY_P],
        bmxpml: vec![vec![0.0; NX_P]; NY_P], bmypml: vec![vec![0.0; NX_P]; NY_P], bmxy: vec![vec![0.0; NX_P]; NY_P], bmyx: vec![vec![0.0; NX_P]; NY_P], bmzx: vec![vec![0.0; NX_P]; NY_P], bmzy: vec![vec![0.0; NX_P]; NY_P],
        expml: vec![vec![0.0; NX_P]; NY_P], eypml: vec![vec![0.0; NX_P]; NY_P], ezx: vec![vec![0.0; NX_P]; NY_P], ezy: vec![vec![0.0; NX_P]; NY_P],
        hxpml: vec![vec![0.0; NX_P]; NY_P], hypml: vec![vec![0.0; NX_P]; NY_P], hzx: vec![vec![0.0; NX_P]; NY_P], hzy: vec![vec![0.0; NX_P]; NY_P], 
        epsd: vec![vec![epsbk; NX_P+1]; NY_P+1], sgmed: vec![vec![sigebk; NX_P+1]; NY_P+1],
        mud: vec![vec![mubk; NX_P+1]; NY_P+1], sgmmd : vec![vec![sigmbk; NX_P+1]; NY_P+1],
        pml_s: Vec::new(), befed: 0.0, duration: 0.0, t0: 0.0, nx: NX, ny: NY};
        //時間ステップ
        let v = C / ((epsbk * mubk).sqrt());
        n.dt = 0.99999/(v * ((1.0 / (DX * DX) + 1.0 / (DY * DY)).sqrt()));

        return n;
    }
}


impl fdtd {
    // 初期化設定(媒質の設定後に実行)
    pub fn setup(&mut self) {
        for y in 0..self.nx as usize{
            for x in 0..self.ny as usize{
                let epsx = 0.5 * (self.epsd[x+1][y+1] + self.epsd[x+1][y]) * EPS0;
                let sgex = 0.5 * (self.sgmed[x+1][y+1] + self.sgmed[x+1][y]);
                let mut a = 0.5 * sgex * self.dt / epsx;
                self.aex[x][y] = (1.0 - a) / (1.0 + a);
                self.bexy[x][y] = self.dt / epsx / (1.0 + a) / DY;

                let epsy = 0.5 * (self.epsd[x+1][y+1] + self.epsd[x][y+1]) * EPS0;
                let sgey = 0.5 * (self.sgmed[x+1][y+1] + self.sgmed[x][y+1]);
                a = 0.5 * sgey * self.dt / epsy;
                self.aey[x][y] = (1.0 - a) / (1.0 + a);
                self.beyx[x][y] = self.dt / epsy / (1.0 + a) / DX;

                let epsz = 0.25 * (self.epsd[x+1][y+1] + self.epsd[x+1][y] + self.epsd[x][y+1] + self.epsd[x][y]) * EPS0;
                let sgez = 0.25 * (self.sgmed[x+1][y+1] + self.sgmed[x+1][y] + self.sgmed[x][y+1] + self.sgmed[x][y]);
                a = 0.5 * sgez * self.dt / epsz;
                self.aez[x][y] = (1.0 - a) / (1.0 + a);
                self.bezy[x][y] = self.dt / epsz / (1.0 + a) / DY;
                self.bezx[x][y] = self.dt / epsz / (1.0 + a) / DX;

                let mux = 0.5 * (self.mud[x+1][y+1] + self.mud[x][y+1]) * MU0;
                let sgmx = 0.5 * (self.sgmmd[x+1][y+1] + self.sgmmd[x][y+1]);
                a = 0.5 * sgmx * self.dt / mux;
                self.amx[x][y] = (1.0 - a) / (1.0 + a);
                self.bmxy[x][y] = self.dt / mux / (1.0 + a) / DY;

                let muy = 0.5 * (self.mud[x+1][y+1] + self.mud[x+1][y]) * MU0;
                let sgmy = 0.5 * (self.sgmmd[x+1][y+1] + self.sgmmd[x+1][y]);
                a = 0.5 * sgmy * self.dt / muy;
                self.amy[x][y] = (1.0 - a) / (1.0 + a);
                self.bmyx[x][y] = self.dt / muy / (1.0 + a) / DX;

                let muz = self.mud[x+1][y+1] * MU0;
                let sgmz = self.sgmmd[x+1][y+1];
                a = 0.5 * sgmz * self.dt / muz;
                self.amz[x][y] = (1.0 - a) / (1.0 + a);
                self.bmzx[x][y] = self.dt / muz / (1.0 + a) / DX;
                self.bmzy[x][y] = self.dt / muz / (1.0 + a) / DY;
            }
        }

        self.init_pml();
    }

    // 電界計算
    pub fn e_cal(&mut self) {
        //Ex
        for y in 1..(self.ny-1) as usize {
            for x in 0..(self.nx-1) as usize {
                self.ex[x][y] = self.aex[x][y] * self.ex[x][y] + self.bexy[x][y] * (self.hz[x][y] - self.hz[x][y-1]);
            }
        }

        //Ey
        for y in 0..(self.ny-1) as usize {
            for x in 1..(self.nx-1) as usize {
                self.ey[x][y] = self.aey[x][y] * self.ey[x][y] - self.beyx[x][y] * (self.hz[x][y] - self.hz[x-1][y]);
            }
        }

        //Ez
        for y in 1..(self.ny-1) as usize {
            for x in 1..(self.nx-1) as usize {
                self.ez[x][y] = self.aez[x][y] * self.ez[x][y] + self.bezx[x][y] * (self.hy[x][y] - self.hy[x-1][y])
                                                               - self.bezy[x][y] * (self.hx[x][y] - self.hx[x][y-1]);
            }
        }
    }

    // 磁界計算
    pub fn h_cal(&mut self) {
        //Hx
        for y in 0..(self.ny-1) as usize {
            for x in 1..(self.nx-1) as usize {
                self.hx[x][y] = self.amx[x][y] * self.hx[x][y] - self.bmxy[x][y] * (self.ez[x][y+1] - self.ez[x][y]);
            }
        }

        //Hy
        for y in 1..(self.ny-1) as usize {
            for x in 0..(self.nx-1) as usize {
                self.hy[x][y] = self.amy[x][y] * self.hy[x][y] + self.bmyx[x][y] * (self.ez[x+1][y] - self.ez[x][y]);
            }
        }

        //Hz
        for y in 0..(self.ny-1) as usize {
            for x in 0..(self.nx-1) as usize {
                self.hz[x][y] = self.amz[x][y] * self.hz[x][y] - self.bmzx[x][y] * (self.ey[x+1][y] - self.ey[x][y])
                                                               + self.bmzy[x][y] * (self.ex[x][y+1] - self.ex[x][y]);
            }
        }
    }

    // PML初期化 1壁
    #[allow(non_snake_case)]
    fn initPml(&mut self, xs: i32, xl: i32, ys: i32, yl: i32) {
        self.pml_s.push(pml{x_s: xs, x_l: xl, y_s: ys, y_l: yl});

        let smax0x = copml * RMAX * (ORDER + 1) as f32 / (LPML as f32 * DX);
        let smax0y = copml * RMAX * (ORDER + 1) as f32 / (LPML as f32 * DY);

        let epspml: f32 = epsbk * EPS0;
        let mupml: f32 = mubk * MU0;

        for y in ys..yl {
            for x in xs..xl {
                let mut sigmxm;
                let mut sigmxe;
                let mut sigmym;
                let mut sigmye;

                let mut a;

                if x < LPML { // 左側のPML初期設定
                    sigmxm = (((LPML - x) as f32 - 0.5) / (LPML as f32)).powi(ORDER) * smax0x;
                    sigmxe = (((LPML - x) as f32) / (LPML as f32)).powi(ORDER) * smax0x;
                }
                else if x >= self.nx - LPML { // 右側のPML初期設定
                    sigmxm = (((x - self.nx + LPML) as f32 + 0.5) / (LPML as f32)).powi(ORDER) * smax0x;
                    sigmxe = (((x - self.nx + LPML) as f32) / (LPML as f32)).powi(ORDER) * smax0x;
                }
                else {
                    sigmxm = 0.0;
                    sigmxe = 0.0;
                }

                if y < LPML { // 上側のPML初期設定
                    sigmym = (((LPML - y) as f32 - 0.5) / (LPML as f32)).powi(ORDER) * smax0y;
                    sigmye = (((LPML - y) as f32) / (LPML as f32)).powi(ORDER) * smax0y;
                }
                else if y >= self.ny - LPML { // 下側のPML初期設定
                    sigmym = (((y - self.ny + LPML) as f32 + 0.5) / (LPML as f32)).powi(ORDER) * smax0y;
                    sigmye = (((y - self.ny + LPML) as f32) / (LPML as f32)).powi(ORDER) * smax0y;
                }
                else {
                    sigmym = 0.0;
                    sigmye = 0.0;
                }

                // PML電界の初期設定
                sigmxe = sigmxe * epsbk;
                a = 0.5 * sigmxe * self.dt / epspml;
                self.aexpml[x as usize][y as usize] = (1.0 - a) / (1.0 + a);
                self.bexpml[x as usize][y as usize] = self.dt / epspml / (1.0 + a) / DX;

                sigmye = sigmye * epsbk;
                a = 0.5 * sigmye * self.dt / epspml;
                self.aeypml[x as usize][y as usize] = (1.0 - a) / (1.0 + a);
                self.beypml[x as usize][y as usize] = self.dt / epspml / (1.0 + a) / DY;

                // PML磁界の初期設定
                sigmxm = sigmxm * epsbk;
                a = 0.5 * sigmxm * self.dt / epspml;
                self.amxpml[x as usize][y as usize] = (1.0 - a) / (1.0 + a);
                self.bmxpml[x as usize][y as usize] = self.dt / mupml / (1.0 + a) / DX;

                sigmym = sigmym * epsbk;
                a = 0.5 * sigmym * self.dt / epspml;
                self.amypml[x as usize][y as usize] = (1.0 - a) / (1.0 + a);
                self.bmypml[x as usize][y as usize] = self.dt / mupml / (1.0 + a) / DY;

            }
        }

    }

    // PML内初期化設定
    fn init_pml(&mut self) {
        let nxt = self.nx.clone();
        let nyt = self.ny.clone();
        self.initPml(0, LPML, 0, nyt.clone());
        self.initPml(nxt.clone()-LPML, nxt.clone(), 0, nyt.clone());
        self.initPml(LPML, nxt.clone()-LPML, 0, LPML);
        self.initPml(LPML, nxt.clone()-LPML, nxt.clone()-LPML, nyt.clone());
    }

    // PML内電界計算
    pub fn e_pml(&mut self){
        let tmp = self.pml_s.clone();
        for n in tmp {
            //Ex
            for y in (n.y_s+1) as usize..(n.y_l-1) as usize {
                for x in (n.x_s as usize)..(n.x_l-1) as usize {
                    self.expml[x][y] = self.aeypml[x][y] * self.expml[x][y] + self.beypml[x][y] * (self.hz[x][y] - self.hz[x][y-1]);
                    self.ex[x][y] = self.expml[x][y].clone();
                }
            }

            //Ey
            for y in (n.y_s as usize)..((n.y_l-1) as usize) {
                for x in ((n.x_s+1) as usize)..((n.x_l-1) as usize) {
                    self.eypml[x][y] = self.aexpml[x][y] * self.eypml[x][y] - self.bexpml[x][y] * (self.hz[x][y] - self.hz[x-1][y]);
                    self.ey[x][y] = self.eypml[x][y].clone();
                }
            }

            //Ez
            for y in ((n.y_s+1) as usize)..((n.y_l-1) as usize) as usize {
                for x in ((n.x_s+1) as usize)..((n.x_l-1) as usize) {
                    self.ezx[x][y] = self.aexpml[x][y] * self.ezx[x][y] + self.bexpml[x][y] * (self.hy[x][y] - self.hy[x-1][y]);
                    self.ezy[x][y] = self.aeypml[x][y] * self.ezy[x][y] - self.beypml[x][y] * (self.hx[x][y] - self.hx[x][y-1]);
                    self.ez[x][y] = self.ezx[x][y] + self.ezy[x][y];
                }
            }
        }

    }

    // PML内磁界計算
    pub fn h_pml(&mut self){
        let tmp = self.pml_s.clone();
        for n in tmp {
            //Hx
            for y in (n.y_s as usize)..((n.y_l - 1) as usize) {
                for x in ((n.x_s + 1) as usize)..((n.x_l - 1) as usize) {
                    self.hxpml[x][y] = self.amypml[x][y] * self.hxpml[x][y] - self.bmypml[x][y] * (self.ez[x][y+1] - self.ez[x][y]);
                    self.hx[x][y] = self.hxpml[x][y].clone();
                }
            }

            //Hy
            for y in ((n.y_s + 1) as usize)..((n.y_l - 1) as usize) {
                for x in (n.x_s as usize)..((n.x_l - 1) as usize) {
                    self.hypml[x][y] = self.amxpml[x][y] * self.hypml[x][y] + self.bmxpml[x][y] * (self.ez[x+1][y] - self.ez[x][y]);
                    self.hy[x][y] = self.hypml[x][y].clone();
                }
            }

            //Hz
            for y in (n.y_s as usize)..((n.y_l - 1) as usize) as usize {
                for x in (n.x_s as usize)..((n.x_l - 1) as usize) {
                    self.hzx[x][y] = self.amxpml[x][y] * self.hzx[x][y] - self.bmxpml[x][y] * (self.ey[x+1][y] - self.ey[x][y]);
                    self.hzy[x][y] = self.amypml[x][y] * self.hzy[x][y] + self.bmypml[x][y] * (self.ex[x][y+1] - self.ex[x][y]);
                    self.hz[x][y] = self.hzx[x][y] + self.hzy[x][y];
                }
            }
        }

    }

    // 障害物媒質設定 epsr: 障害物の誘電率
    pub fn epsmu(&mut self, x_s: usize, x_l: usize, y_s: usize, y_l: usize, epsr: f32){
        for y in y_s+1..y_l {
            for x in x_s+1..x_l {
                self.epsd[x][y] = epsr;
                self.mud[x][y] = 1.0;
                self.sgmed[x][y] = 0.0;
                self.sgmmd[x][y] = 0.0;
            }
        }
    }

    // 障害物設定 完全導体
    pub fn pec_rect(&mut self, x_s: usize, x_l: usize, y_s: usize, y_l: usize){
        for y in y_s..y_l {
            for x in x_s..x_l-1 {
                self.aex[x][y] = 0.0;
                self.bexy[x][y] = 0.0;
            }
        }
        for y in y_s..y_l-1 {
            for x in x_s..x_l {
                self.aey[x][y] = 0.0;
                self.beyx[x][y] = 0.0;
            }
        }
        for y in y_s..y_l {
            for x in x_s..x_l {
                self.aez[x][y] = 0.0;
                self.bezx[x][y] = 0.0;
                self.bezy[x][y] = 0.0;
            }
        }
    }

    // 電流源の初期化
    pub fn init_source(&mut self, x: usize, y: usize){
        let epsz = 0.25 * (self.epsd[x+1][y+1] + self.epsd[x][y+1] + self.epsd[x+1][y] + self.epsd[x][y]) * EPS0;
        self.befed = self.dt / epsz;
        self.duration = 0.1e-9;
        self.t0 = self.duration * 4.0;
    }

    // 電流の計算
    pub fn feed(&mut self, x: usize, y: usize, t: f32){
        let mut tmp:f32 = (t - 0.5 * self.dt - self.t0) / self.duration;
        tmp = tmp.powi(2);
        let iz = (-tmp).exp();
        //let iz = (-((t - 0.5 * self.dt - self.t0) / self.duration).powf(2.0)).exp();
        self.ez[x][y] = self.ez[x][y] - self.befed * iz / (DX * DY);
        println!("ez:{} iz:{} befed:{}", self.ez[x][y], iz, self.befed);
    }

    pub fn out_p(&self, xo: usize, yo: usize) {
        for y in LPML as usize..(self.ny - LPML) as usize {
            for x in LPML as usize.. (self.nx - LPML) as usize {
                println!("ez[{}][{}] = {}", x-LPML as usize, y-LPML as usize, self.ez[x][y]);
            }
        }
        println!("Observation point: {}",self.ez[xo][yo]);
    }

    #[allow(unused_must_use)]
    pub fn out_file(&self, file: &mut File, xo: usize, yo: usize) {
        for y in LPML as usize..(self.ny - LPML) as usize {
            for x in LPML as usize.. (self.nx - LPML) as usize {
                write!(*file, "ez[{}][{}] = {}\r\n", x-LPML as usize, y-LPML as usize, self.ez[x][y]);
            }
        }
        write!(*file,"Observation point: {}\r\n",self.ez[xo][yo]);
    }
}
