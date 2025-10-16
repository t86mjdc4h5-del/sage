# -*- coding: utf-8 -*-
# bsd_enhanced_collector.sage - النسخة المحسنة جاهزة للتشغيل
from sage.all import *
from sage.parallel.decorate import parallel
import csv, time, os
from sage.databases.cremona import CremonaDatabase

class BSDDataCollector:
    def __init__(self, precision=80):
        self.precision = precision
        self.results = []
        
    def get_curve_data(self, label):
        """Collect comprehensive BSD data for a single curve"""
        try:
            print(f"Processing {label}...")
            E = EllipticCurve(label)
            
            # Basic invariants
            N = E.conductor()
            r = E.rank()
            tors = E.torsion_subgroup().cardinality()
            
            # Precision setup
            RR = RealField(self.precision)
            
            # Period computation with component awareness
            lattice = E.period_lattice()
            try:
                # lattice.omega() may return one period; adjust for real components
                if E.real_components() == 2:
                    Omega = RR(2 * lattice.omega())
                else:
                    Omega = RR(lattice.omega())
            except Exception:
                # Fallback to E.real_periods()
                try:
                    real_periods = E.real_periods()
                    Omega = RR(real_periods[0]) if real_periods else RR(0)
                except Exception:
                    Omega = RR(0)
            
            # Regulator - handle rank 0 case
            try:
                Reg = RR(E.regulator()) if r > 0 else RR(1)
            except Exception:
                Reg = RR(1) if r == 0 else RR(0)
            
            # Tamagawa product
            tam_prod = 1
            try:
                for p in E.bad_primes():
                    tam_prod *= E.tamagawa_number(p)
            except Exception:
                tam_prod = 1
            
            # L-series computation with precision (dokchitser backend)
            try:
                L = E.lseries().dokchitser(prec=self.precision)
            except Exception:
                # fallback to LSeries without dokchitser if unavailable
                L = E.lseries()
            
            if r == 0:
                try:
                    L_val = RR(L(1))
                except Exception:
                    L_val = RR(0)
                L_derivative = RR(0)
            else:
                # Compute r-th derivative at s=1; use L.derivative(k)(1) if supported
                try:
                    # Some lseries implementations: L.derivative(order)(s)
                    Lr = L.derivative(r)(1)
                    L_derivative = RR(Lr / factorial(r))
                    L_val = RR(Lr)
                except Exception:
                    # fallback: try different API
                    try:
                        Lr = L.derivative(1, r)
                        L_derivative = RR(Lr / factorial(r))
                        L_val = RR(Lr)
                    except Exception:
                        L_val = RR(0)
                        L_derivative = RR(0)
            
            # Analytic rank from L-series if available
            analytic_rank = None
            try:
                analytic_rank = int(L.rank())
            except Exception:
                analytic_rank = None
            
            # BSD ratio calculation
            denominator = Omega * Reg * tam_prod / (tors**2) if tors != 0 else Omega * Reg * tam_prod
            try:
                if denominator != 0:
                    if r > 0:
                        bsd_ratio = RR(L_derivative / denominator)
                    else:
                        bsd_ratio = RR(L_val / denominator)
                else:
                    bsd_ratio = RR(0)
            except Exception:
                bsd_ratio = RR(0)
            
            # Heegner index if applicable (for rank 1 curves)
            heegner_index = None
            if r == 1 and E.root_number() == -1:
                try:
                    # This might take time for some curves
                    P = E.heegner_point(-1, prec=max(80, self.precision))
                    if P is not None and hasattr(P, 'is_zero') and not P.is_zero():
                        h_pt = P.height()
                        h_reg = Reg
                        if h_reg > 0:
                            heegner_index = RR(h_pt / h_reg)
                except Exception:
                    heegner_index = None
            
            return {
                'label': label,
                'conductor': N,
                'rank': r,
                'analytic_rank': analytic_rank,
                'root_number': E.root_number(),
                'L_value': float(L_val) if L_val is not None else None,
                'L_derivative': float(L_derivative) if L_derivative is not None else None,
                'Omega': float(Omega) if Omega is not None else None,
                'Regulator': float(Reg) if Reg is not None else None,
                'Tamagawa_product': tam_prod,
                'torsion_size': tors,
                'bsd_ratio': float(bsd_ratio) if bsd_ratio is not None else None,
                'heegner_index': float(heegner_index) if heegner_index else None,
                'real_components': E.real_components(),
                'cm_discriminant': E.cm_discriminant() if E.has_cm() else 0,
                'status': 'success'
            }
            
        except Exception as e:
            print(f"Error processing {label}: {str(e)}")
            return {'label': label, 'status': f'error: {str(e)}'}
    
    def collect_curated_sample(self, sample_labels):
        """Collect data for a curated sample of curves"""
        print(f"🎯 Starting curated sample collection for {len(sample_labels)} curves...")
        
        for i, label in enumerate(sample_labels, 1):
            print(f"Progress: {i}/{len(sample_labels)} — {label}")
            data = self.get_curve_data(label)
            self.results.append(data)
    
    def collect_sample(self, sample_size=20, max_conductor=10000):
        """Collect data for a systematic sample"""
        c = CremonaDatabase()
        curves_collected = 0
        
        # We'll sample across different conductor ranges and ranks
        conductor_ranges = [(1, 1000), (1001, 5000)]
        target_ranks = [0, 1, 2]  # Focus on these ranks initially
        
        for cond_range in conductor_ranges:
            if curves_collected >= sample_size:
                break
                
            for N in range(cond_range[0], cond_range[1] + 1):
                if curves_collected >= sample_size:
                    break
                    
                try:
                    curves = c.curves(N)
                    for curve_info in curves:
                        if curves_collected >= sample_size:
                            break
                            
                        label = curve_info['label']
                        data = self.get_curve_data(label)
                        
                        # Only keep curves with ranks in our target
                        if data.get('rank') in target_ranks:
                            self.results.append(data)
                            curves_collected += 1
                            print(f"Collected {curves_collected}/{sample_size}: {label}")
                            
                except Exception:
                    continue
    
    def save_to_csv(self, filename='bsd_enhanced_data.csv'):
        """Save results to CSV"""
        if not self.results:
            print("No data to save!")
            return
            
        fieldnames = ['label', 'conductor', 'rank', 'analytic_rank', 'root_number', 
                     'L_value', 'L_derivative', 'Omega', 'Regulator', 'Tamagawa_product',
                     'torsion_size', 'bsd_ratio', 'heegner_index', 'real_components', 
                     'cm_discriminant', 'status']
        
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            for result in self.results:
                writer.writerow(result)
        
        print(f"Data saved to {filename}")

# العينة المختارة بعناية للاختبار الأولي
CURATED_SAMPLE = [
    # Rank 0 cases
    "11a1", "37a1", "43a1", "53a1", "61a1", 
    
    # Rank 1 cases  
    "37a1", "389a1", "5077a1", "110a1", "210e1",
    
    # Rank 2 cases
    "389a1", "431a1", "446d1", "571a1", "681b1",
    
    # Interesting cases
    "571a1", "681b1", "2773a1", "446d1",
    
    # CM curves
    "27a1", "32a1", "36a1", "49a1"
]

# التنفيذ الفوري
if __name__ == "__main__":
    # إنشاء المجلد المخصص
    os.makedirs('bsd_data', exist_ok=True)
    os.chdir('bsd_data')
    
    collector = BSDDataCollector(precision=80)
    
    print("🚀 Starting enhanced BSD data collection...")
    print("📁 Working directory:", os.getcwd())
    
    start_time = time.time()
    
    # المرحلة 1: تشغيل العينة المختارة (20 منحنى)
    collector.collect_curated_sample(CURATED_SAMPLE)
    collector.save_to_csv('bsd_enhanced_data_phase1.csv')
    
    elapsed = time.time() - start_time
    print(f"✅ Phase 1 completed in {elapsed:.2f} seconds")
    print(f"📊 Collected data for {len(collector.results)} curves")
    
    # تحليل أولي فوري
    print("\n📈 IMMEDIATE ANALYSIS:")
    success_data = [r for r in collector.results if r.get('status') == 'success']
    print(f"Successful collections: {len(success_data)}")
    
    for data in success_data[:5]:  # عرض أول 5 نتائج للتحقق
        try:
            print(f"{data['label']}: rank={data['rank']}, L_value={data['L_value']:.6f}, "
                  f"Omega={data['Omega']:.6f}, bsd_ratio={data['bsd_ratio']:.6f}")
        except Exception:
            print(f"{data['label']}: (some values missing) status={data.get('status')}")