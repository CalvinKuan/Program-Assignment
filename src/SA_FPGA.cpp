#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <math.h>
#include <random>
#include <chrono>
#include <unordered_set>
#include <algorithm>
#include <climits>

std::mt19937 rng(123456);

using namespace std;
using ns = std::chrono::nanoseconds;
using Clock = std::chrono::high_resolution_clock;

int R = 0, C = 0;

struct pin_info
{
    string name;
    double x;
    double y;
};

struct logic_info
{
    string name;
    int x;
    int y;
    vector<int> connect_nets;
};

struct rec
{
    double xmax;
    double xmin;
    double ymax;
    double ymin;
};

// net 端點：是 block 還是 pin
struct NetTerm
{
    bool is_pin; // false => logic block, true => IO pin
    int idx;     // index in logic_blocks_table 或 pin_table
};

struct net_info
{
    string name;
    int degree;
    vector<string> connections;  // 原始名稱（保留給 debug）
    vector<NetTerm> terms;       // 給 HPWL 用，不再查 hash
};

struct net_hpwl_cong
{
    double hpwl;
    rec box;
};

struct change_info
{
    string nameA, nameB;
    int Ax, Ay;
    int Bx, By;
};

// 給 incremental 用的小結構
struct NetUpdate
{
    int net_idx;
    std::string net_name;
    double old_hpwl;
    rec old_box;
    double new_hpwl;
    rec new_box;
};

struct CellBackup
{
    int x, y;
    int old_val;
};

void initial_place(vector<vector<string>> &array2D, int R, int C,
                   vector<logic_info> &logic_blocks_table)
{
    int N = (int)logic_blocks_table.size();

    int idx = 0;
    for (int i = 0; i < R && idx < N; ++i)
    {
        for (int j = 0; j < C && idx < N; ++j)
        {
            array2D[i][j] = logic_blocks_table[idx].name;
            logic_blocks_table[idx].x = i;  // x = row
            logic_blocks_table[idx].y = j;  // y = col
            ++idx;
        }
    }
}

// 新版：不再查 name2idx/pin2idx，只看 net.terms
rec calculate_bounding_box(const net_info &net,
                           const vector<logic_info> &logic_blocks_table,
                           const vector<pin_info> &pin_table)
{
    bool first = true;
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;

    for (const auto &term : net.terms)
    {
        if (term.idx < 0) continue;

        double this_xmin, this_xmax, this_ymin, this_ymax;

        if (!term.is_pin)
        {
            // movable block
            const auto &blk = logic_blocks_table[term.idx];
            double xb = (double)blk.x;
            double yb = (double)blk.y;
            this_xmin = xb;
            this_xmax = xb + 1.0;
            this_ymin = yb;
            this_ymax = yb + 1.0;
        }
        else
        {
            // fixed pin
            const auto &p = pin_table[term.idx];
            double xp = p.x;
            double yp = p.y;
            this_xmin = xp;
            this_xmax = xp;
            this_ymin = yp;
            this_ymax = yp;
        }

        if (first)
        {
            xmin = this_xmin;
            xmax = this_xmax;
            ymin = this_ymin;
            ymax = this_ymax;
            first = false;
        }
        else
        {
            xmin = min(xmin, this_xmin);
            xmax = max(xmax, this_xmax);
            ymin = min(ymin, this_ymin);
            ymax = max(ymax, this_ymax);
        }
    }

    rec Rect;
    if (first)
    {
        Rect.xmin = Rect.xmax = Rect.ymin = Rect.ymax = 0.0;
        return Rect;
    }

    Rect.xmin = xmin;
    Rect.xmax = xmax;
    Rect.ymin = ymin;
    Rect.ymax = ymax;
    return Rect;
}

void compute_hpwl_cong(const vector<logic_info> &logic_blocks_table,
                       const vector<pin_info> &pin_table,
                       const vector<net_info> &net_table,
                       vector<vector<int>> &U,
                       unordered_map<string, net_hpwl_cong> &net2hpwl_cong,
                       double &total_hpwl,
                       double &CC)
{
    // 初始化 coverage
    for (int y = 0; y < R; ++y)
    {
        fill(U[y].begin(), U[y].end(), 0);
    }

    total_hpwl = 0.0;

    for (size_t i = 0; i < net_table.size(); ++i)
    {
        const net_info &net = net_table[i];
        rec box = calculate_bounding_box(net, logic_blocks_table, pin_table);
        double hpwl = (box.xmax - box.xmin) + (box.ymax - box.ymin);
        total_hpwl += hpwl;

        net_hpwl_cong &info = net2hpwl_cong[net.name];
        info.hpwl = hpwl;
        info.box = box;

        int x_start = max(0, (int)floor(box.xmin));
        int x_end   = min(C, (int)ceil(box.xmax)); // [x_start, x_end)
        int y_start = max(0, (int)floor(box.ymin));
        int y_end   = min(R, (int)ceil(box.ymax)); // [y_start, y_end)

        for (int y = y_start; y < y_end; ++y)
        {
            for (int x = x_start; x < x_end; ++x)
            {
                U[y][x] += 1;
            }
        }
    }

    // 計算 CC
    long long N = (long long)R * (long long)C;
    long long sumU = 0;
    long long sumU2 = 0;
    for (int y = 0; y < R; ++y)
    {
        for (int x = 0; x < C; ++x)
        {
            int v = U[y][x];
            sumU += v;
            sumU2 += 1LL * v * v;
        }
    }

    if (sumU == 0)
    {
        CC = 1.0; // 沒有任何 net 覆蓋時，視為完美均勻
    }
    else
    {
        CC = (double)N * (double)sumU2 / ((double)sumU * (double)sumU);
    }
}

double compute_cost(double total_hpwl, double CC, double /*lambda*/)
{
    // 目前直接用 total_hpwl * CC
    return total_hpwl * CC;
}

void undo_move(vector<logic_info> &logic_blocks_table,
               unordered_map<string, int> &name2idx,
               vector<vector<string>> &array2D,
               const change_info &chg)
{
    // 還原 A
    int idxA = name2idx[chg.nameA];
    auto &A = logic_blocks_table[idxA];

    if (chg.nameB != "none")
    {
        // block <-> block 的情況
        int idxB = name2idx[chg.nameB];
        auto &B = logic_blocks_table[idxB];

        A.x = chg.Ax;
        A.y = chg.Ay;
        B.x = chg.Bx;
        B.y = chg.By;

        array2D[chg.Ax][chg.Ay] = chg.nameA;
        array2D[chg.Bx][chg.By] = chg.nameB;
    }
    else
    {
        // block <-> 空格 的情況
        A.x = chg.Ax;
        A.y = chg.Ay;

        array2D[chg.Ax][chg.Ay] = chg.nameA;
        array2D[chg.Bx][chg.By] = "none";
    }
}

// ========= Optimal Region 版 random move =========
change_info random_move_optimal(vector<logic_info> &logic_blocks_table,
                                unordered_map<string, int> &name2idx,
                                vector<vector<string>> &array2D,
                                int R, int C,
                                const vector<net_info> &net_table,
                                const unordered_map<string, net_hpwl_cong> &net2hpwl_cong)
{
    int n = (int)logic_blocks_table.size();
    std::uniform_int_distribution<int> dist_block(0, n - 1);
    change_info changes;

    // 1. 隨機選一個 block A
    int idx = dist_block(rng);
    auto &A = logic_blocks_table[idx];
    int Ax = A.x;
    int Ay = A.y;

    int Bx, By;

    // 2. 依據 A 連到的 nets 算「目標區域」中心點
    double cx_sum = 0.0, cy_sum = 0.0;
    int cnt = 0;

    for (int nid : A.connect_nets)
    {
        const net_info &net = net_table[nid];
        auto it_info = net2hpwl_cong.find(net.name);
        if (it_info == net2hpwl_cong.end()) continue;

        const rec &b = it_info->second.box;
        double cx = 0.5 * (b.xmin + b.xmax);
        double cy = 0.5 * (b.ymin + b.ymax);
        cx_sum += cx;
        cy_sum += cy;
        cnt++;
    }

    bool use_region = (cnt > 0);

    if (use_region)
    {
        double cx = cx_sum / (double)cnt;
        double cy = cy_sum / (double)cnt;

        // 在 net 中心附近開一個 window 當 optimal region
        int half_range_r = std::max(2, R / 10); // row 方向
        int half_range_c = std::max(2, C / 10); // col 方向

        int r_min = std::max(0, (int)floor(cx - half_range_r));
        int r_max = std::min(R - 1, (int)ceil (cx + half_range_r));
        int c_min = std::max(0, (int)floor(cy - half_range_c));
        int c_max = std::min(C - 1, (int)ceil (cy + half_range_c));

        if (r_min > r_max || c_min > c_max)
        {
            use_region = false; // 算壞掉就用全域隨機
        }
        else
        {
            std::uniform_int_distribution<int> distR(r_min, r_max);
            std::uniform_int_distribution<int> distC(c_min, c_max);

            Bx = distR(rng);
            By = distC(rng);

            // 萬一 region 只有自己那一格，就 fallback
            if (Bx == Ax && By == Ay)
            {
                use_region = false;
            }
        }
    }

    // 3. 如果沒辦法用 optimal region，就回到原本「整張板子亂丟」
    if (!use_region)
    {
        std::uniform_int_distribution<int> distR(0, R - 1);
        std::uniform_int_distribution<int> distC(0, C - 1);
        Bx = distR(rng);
        By = distC(rng);
        while (Bx == Ax && By == Ay)
        {
            Bx = distR(rng);
            By = distC(rng);
        }
    }

    // 4. 做 swap / move
    string swap_target = array2D[Bx][By];
    string temp;

    if (swap_target != "none")
    {
        // block <-> block
        int idxB = name2idx.at(swap_target);
        auto &B = logic_blocks_table[idxB];
        temp = B.name;

        A.x = Bx; A.y = By;
        B.x = Ax; B.y = Ay;

        array2D[Ax][Ay] = swap_target;
        array2D[Bx][By] = A.name;
    }
    else
    {
        // block <-> 空格
        A.x = Bx;
        A.y = By;
        temp = "none";

        array2D[Bx][By] = A.name;
        array2D[Ax][Ay] = "none";
    }

    changes.nameA = A.name;
    changes.Ax = Ax;
    changes.Ay = Ay;
    changes.nameB = temp;
    changes.Bx = Bx;
    changes.By = By;
    return changes;
}

// =================================================

void SA(vector<logic_info> &logic_blocks_table,
        vector<vector<string>> &array2D,
        vector<pin_info> &pin_table,
        vector<net_info> &net_table,
        unordered_map<string, int> &name2idx,
        unordered_map<string, int> &pin2idx, // 保留參數（雖然目前 bbox 不再用）
        unordered_map<string, net_hpwl_cong> &net2hpwl_cong,
        vector<vector<int>> &U,
        double &total_hpwl,
        double &CC)
{
    auto start_SA = Clock::now();
    const int TIME_LIMIT_MS = 220000; // 仍然保留 220 秒安全網

    int num_cells = (int)logic_blocks_table.size();
    int num_nets  = (int)net_table.size();

    double pen0 = std::max(0.0, CC - 1.0);
    if (pen0 < 1e-6)
        pen0 = 1.0; // 避免除 0
    double lambda = 0.15 * total_hpwl / pen0;
    double cur_cost = compute_cost(total_hpwl, CC, lambda);
    double best_cost = cur_cost;
    // 紀錄最佳解
    auto best_logic_blocks = logic_blocks_table;
    auto best_array2D = array2D;
    double best_hpwl = total_hpwl;
    double best_CC = CC;

    cerr << "[INIT] HPWL=" << total_hpwl
         << " CC=" << CC
         << " cost=" << cur_cost << endl;

    // ===== SA 參數 (調過版) =====
    double T      = 500.0;
    double T_end  = 0.0001;
    double alpha  = 0.95;

    // 根據 problem size 自動調迭代次數
    int base_iter = std::max(2000, std::min(num_cells * 2, 20000));

    std::uniform_real_distribution<double> dist01(0.0, 1.0);

    // 先從現在的 U 算一次 sumU / sumU2，之後都 incremental 更新
    const long long N = (long long)R * (long long)C;
    long long sumU = 0;
    long long sumU2 = 0;
    for (int y = 0; y < R; ++y)
    {
        for (int x = 0; x < C; ++x)
        {
            int v = U[y][x];
            sumU += v;
            sumU2 += 1LL * v * v;
        }
    }

    // 使用一維 visit_tag 陣列取代 unordered_map visited（U 還原用）
    const int cellN = R * C;
    std::vector<int> visit_tag(cellN, 0);
    int cur_tag = 0;

    // net_mark 陣列取代 unordered_set affected_net_set
    int netN = (int)net_table.size();
    std::vector<int> net_mark(netN, 0);
    int net_tag = 0;

    // 早停控制：連續幾個溫度沒有明顯改善就停
    int    stall_temps     = 0;
    double last_best_cost  = best_cost;
    const  double IMPROVE_EPS = 1e-3;  // 視為有改善的 threshold
    const  int    STALL_LIMIT = 20;    // 連續 20 個溫度沒進步就停

    while (T > T_end)
    {
        auto now = Clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - start_SA);
        if (duration.count() > TIME_LIMIT_MS)
        {
            cerr << "[SA] Time limit reached, break.\n";
            break;
        }

        if (stall_temps >= STALL_LIMIT)
        {
            cerr << "[SA] Early stop: no improvement for " << STALL_LIMIT << " temperature steps.\n";
            break;
        }

        int iter_per_T = base_iter;

        for (int it = 0; it < iter_per_T; ++it)
        {
            // 更新 tag
            cur_tag++;
            if (cur_tag == INT_MAX)
            {
                std::fill(visit_tag.begin(), visit_tag.end(), 0);
                cur_tag = 1;
            }

            net_tag++;
            if (net_tag == INT_MAX)
            {
                std::fill(net_mark.begin(), net_mark.end(), 0);
                net_tag = 1;
            }

            // 1. 做一個 optimal region random move
            change_info chg = random_move_optimal(
                logic_blocks_table, name2idx,
                array2D, R, C,
                net_table, net2hpwl_cong
            );

            // 2. 找出受影響的 nets（A、B 兩顆 cell 相關聯的 nets）
            std::vector<int> affected_nets;
            affected_nets.reserve(32);

            auto collect_nets = [&](const std::string &blk_name)
            {
                if (blk_name == "none")
                    return;
                auto it_blk = name2idx.find(blk_name);
                if (it_blk == name2idx.end())
                    return;
                const auto &blk = logic_blocks_table[it_blk->second];
                for (int nid : blk.connect_nets)
                {
                    if (net_mark[nid] != net_tag)
                    {
                        net_mark[nid] = net_tag;
                        affected_nets.push_back(nid);
                    }
                }
            };

            collect_nets(chg.nameA);
            collect_nets(chg.nameB);

            // 3. 只對 affected nets 用 net2hpwl_cong 查舊 hpwl / box，
            //    然後重新算新的 bbox / hpwl
            std::vector<NetUpdate> updates;
            updates.reserve(affected_nets.size());

            double trial_hpwl = total_hpwl; // 從現在的 total_hpwl 開始加差值

            for (int nid : affected_nets)
            {
                const net_info &net = net_table[nid];
                const std::string &net_name = net.name;

                // 從 net2hpwl_cong 查舊 hpwl / box
                auto it_info = net2hpwl_cong.find(net_name);
                double old_hp = 0.0;
                rec old_rect{};
                if (it_info != net2hpwl_cong.end())
                {
                    old_hp = it_info->second.hpwl;
                    old_rect = it_info->second.box;
                }

                // 用目前 (已 move) 的 logic_blocks_table/pin_table
                // 重新算 bbox
                rec new_rect = calculate_bounding_box(net,
                                                      logic_blocks_table,
                                                      pin_table);
                double new_hp = (new_rect.xmax - new_rect.xmin) +
                                (new_rect.ymax - new_rect.ymin);

                trial_hpwl += (new_hp - old_hp);

                NetUpdate nu;
                nu.net_idx = nid;
                nu.net_name = net_name;
                nu.old_hpwl = old_hp;
                nu.old_box = old_rect;
                nu.new_hpwl = new_hp;
                nu.new_box = new_rect;
                updates.push_back(nu);
            }

            // 4. 對 U 做「暫時性的」更新，順便算 new_sumU / new_sumU2
            long long new_sumU = sumU;
            long long new_sumU2 = sumU2;

            std::vector<CellBackup> cell_backups;
            cell_backups.reserve(256);

            auto apply_rect_delta = [&](const rec &box, int delta)
            {
                int x_start = std::max(0, (int)floor(box.xmin));
                int x_end   = std::min(C, (int)ceil(box.xmax)); // [x_start, x_end)
                int y_start = std::max(0, (int)floor(box.ymin));
                int y_end   = std::min(R, (int)ceil(box.ymax)); // [y_start, y_end)

                for (int y = y_start; y < y_end; ++y)
                {
                    for (int x = x_start; x < x_end; ++x)
                    {
                        int id = y * C + x;

                        if (visit_tag[id] != cur_tag)
                        {
                            visit_tag[id] = cur_tag;
                            CellBackup bk{x, y, U[y][x]};
                            cell_backups.push_back(bk);
                        }

                        int old_v = U[y][x];
                        int new_v = old_v + delta;
                        U[y][x] = new_v;

                        new_sumU  += delta;
                        new_sumU2 += 1LL * new_v * new_v - 1LL * old_v * old_v;
                    }
                }
            };

            // 用 net2hpwl_cong 的舊 box 先 -1，再用新 box +1
            for (const auto &u : updates)
            {
                if (u.old_hpwl > 0.0)
                {
                    apply_rect_delta(u.old_box, -1);
                }
                apply_rect_delta(u.new_box, +1);
            }

            // 用更新後的 U 統計值得到 new_CC
            double new_CC;
            if (new_sumU == 0)
            {
                new_CC = 1.0;
            }
            else
            {
                new_CC = (double)N * (double)new_sumU2 /
                         ((double)new_sumU * (double)new_sumU);
            }

            double new_cost = compute_cost(trial_hpwl, new_CC, lambda);
            double delta_cost = new_cost - cur_cost;
            bool accept = false;

            // 5. SA 接受規則
            if (delta_cost <= 0)
            {
                accept = true;
            }
            else
            {
                double u = dist01(rng);
                if (u < exp(-delta_cost / T))
                    accept = true;
            }

            if (accept)
            {
                // 6. 接受：commit incremental 更新
                cur_cost = new_cost;
                total_hpwl = trial_hpwl;
                CC = new_CC;
                sumU = new_sumU;
                sumU2 = new_sumU2;

                // 用 net2hpwl_cong 更新 affected nets 的 hpwl / box
                for (const auto &u : updates)
                {
                    auto &info = net2hpwl_cong[u.net_name];
                    info.hpwl = u.new_hpwl;
                    info.box = u.new_box;
                }

                // 更新 best solution
                if (cur_cost < best_cost)
                {
                    best_cost = cur_cost;
                    best_hpwl = total_hpwl;
                    best_CC = CC;
                    best_logic_blocks = logic_blocks_table;
                    best_array2D = array2D;
                }
            }
            else
            {
                // 7. 不接受：把 U 還原，位置用 undo_move 拉回
                for (const auto &bk : cell_backups)
                {
                    U[bk.y][bk.x] = bk.old_val;
                }
                undo_move(logic_blocks_table, name2idx, array2D, chg);
                // total_hpwl / CC / sumU / sumU2 / net2hpwl_cong 都維持舊的
            }
        }

        cerr << "[TEMP] T=" << T
             << " best_cost=" << best_cost
             << " best_HPWL=" << best_hpwl
             << " best_CC=" << best_CC << endl;

        // 溫度層結束後檢查是否有進步
        if (best_cost < last_best_cost - IMPROVE_EPS)
        {
            last_best_cost = best_cost;
            stall_temps = 0;
        }
        else
        {
            stall_temps++;
        }

        T *= alpha; // 降溫
    }

    // 退火結束 → 回到最佳解
    logic_blocks_table = best_logic_blocks;
    array2D = best_array2D;
    total_hpwl = best_hpwl;
    CC = best_CC;

    // 最後再 full recompute 一次，讓 net2hpwl_cong / U 對應到 best solution
    compute_hpwl_cong(logic_blocks_table, pin_table, net_table,
                      U, net2hpwl_cong, total_hpwl, CC);
    best_cost = compute_cost(total_hpwl, CC, lambda);

    cerr << "[FINAL] HPWL=" << total_hpwl
         << " CC=" << CC
         << " cost=" << best_cost << endl;
}


void parse_input(ifstream &parse, vector<logic_info> &logic_blocks_table, vector<pin_info> &pin_table, vector<net_info> &net_table, unordered_map<string, int> &name2idx, unordered_map<string, int> &pin2idx, unordered_map<string, int> &net2idx)
{
    string line;
    getline(parse, line);
    istringstream info(line);
    int logic_block_R, logic_block_C, num_logic_blocks, num_io_pins, num_nets;
    info >> logic_block_R >> logic_block_C >> num_logic_blocks >> num_io_pins >> num_nets;
    R = logic_block_R;
    C = logic_block_C;
    for (int i = 0; i < num_logic_blocks; i++)
    {
        string name;
        getline(parse, line);
        istringstream logic_inf(line);
        logic_inf >> name;
        logic_info inst;
        inst.name = name;
        inst.x = 0;
        inst.y = 0;
        logic_blocks_table.push_back(inst);
        name2idx.insert({name, i});
    }
    for (int i = 0; i < num_io_pins; i++)
    {
        string name;
        double x, y;
        getline(parse, line);
        istringstream pin_inf(line);
        pin_inf >> name >> x >> y;
        pin_info inst;
        inst.name = name;
        inst.x = x;
        inst.y = y;
        pin_table.push_back(inst);
        pin2idx.insert({name, i});
    }
    for (int i = 0; i < num_nets; i++)
    {
        string name;
        int degree;
        getline(parse, line);
        istringstream net_inf(line);
        net_inf >> name >> degree;
        net_info inst;
        inst.name = name;
        inst.degree = degree;

        for (int j = 0; j < degree; j++)
        {
            string item;
            net_inf >> item;
            inst.connections.push_back(item);

            NetTerm term;
            auto it_blk = name2idx.find(item);
            if (it_blk != name2idx.end())
            {
                // 是 movable block
                term.is_pin = false;
                term.idx = it_blk->second;
                logic_blocks_table[it_blk->second].connect_nets.push_back(i);
            }
            else
            {
                auto it_pin = pin2idx.find(item);
                if (it_pin != pin2idx.end())
                {
                    // 是 IO pin
                    term.is_pin = true;
                    term.idx = it_pin->second;
                }
                else
                {
                    // 理論上不會
                    term.is_pin = true;
                    term.idx = -1;
                }
            }
            inst.terms.push_back(term);
        }
        net_table.push_back(inst);
        net2idx.insert({name, i});
    }
}

void print_logic_blocks(const vector<logic_info> &logic_blocks_table)
{
    cout << "=== Logic Blocks ===\n";
    for (size_t i = 0; i < logic_blocks_table.size(); ++i)
    {
        const auto &lb = logic_blocks_table[i];
        cout << "[" << i << "] "
             << "name = " << lb.name
             << ", pos = (" << lb.x << ", " << lb.y << ")"
             << ", nets = { ";

        for (size_t k = 0; k < lb.connect_nets.size(); ++k)
        {
            cout << lb.connect_nets[k];
            if (k + 1 != lb.connect_nets.size())
                cout << ", ";
        }

        cout << " }\n";
    }
    cout << endl;
}
void print_pins(const vector<pin_info> &pin_table)
{
    cout << "=== IO Pins ===\n";
    for (size_t i = 0; i < pin_table.size(); ++i)
    {
        const auto &p = pin_table[i];
        cout << i << ": name=" << p.name
             << " x=" << p.x
             << " y=" << p.y << "\n";
    }
    cout << endl;
}
void print_nets(const vector<net_info> &net_table)
{
    cout << "=== Nets ===\n";
    for (size_t i = 0; i < net_table.size(); ++i)
    {
        const auto &n = net_table[i];
        cout << i << ": name=" << n.name
             << " degree=" << n.degree
             << " connections = ";
        for (auto &c : n.connections)
        {
            cout << c << " ";
        }
        cout << "\n";
    }
    cout << endl;
}
void print_net2hpwl_cong(const unordered_map<string, net_hpwl_cong> &net2hpwl_cong)
{
    cout << "=== net2hpwl_cong ===\n";
    int count = 0;
    for (const auto &kv : net2hpwl_cong)
    {
        const string &net_name = kv.first;
        const net_hpwl_cong &info = kv.second;

        cout << "Net: " << net_name << "\n";
        cout << "  HPWL = " << info.hpwl << "\n";
        cout << "  Box: "
             << "(xmin=" << info.box.xmin
             << ", xmax=" << info.box.xmax
             << ", ymin=" << info.box.ymin
             << ", ymax=" << info.box.ymax
             << ")\n";
        count++;
    }
    cout << count << endl;
    cout << endl;
}
void print_U(const vector<vector<int>> &U)
{
    cout << "=== U congestion map ===\n";
    for (int r = 0; r < (int)U.size(); r++)
    {
        for (int c = 0; c < (int)U[r].size(); c++)
        {
            cout << U[r][c] << " ";
        }
        cout << "\n";
    }
    cout << endl;
}
void print_array2D(const vector<vector<string>> &array2D)
{
    cout << "=== array2D placement ===\n";
    for (int r = 0; r < (int)array2D.size(); r++)
    {
        for (int c = 0; c < (int)array2D[r].size(); c++)
        {
            cout << array2D[r][c] << " ";
        }
        cout << "\n";
    }
    cout << endl;
}
void write_output(const string &out_file, const vector<logic_info> &logic_blocks_table)
{
    ofstream fout(out_file);
    if (!fout)
    {
        cerr << "Failed to open output file: " << out_file << endl;
        return;
    }
    for (const auto &lb : logic_blocks_table)
    {
        fout << lb.name << " " << lb.x << " " << lb.y << "\n";
    }
    fout.close();
}

int main(int argc, char *argv[])
{
    auto start = Clock::now();

    if (argc < 3)
    {
        cerr << "Usage: " << argv[0] << " <input_file> <out_file>" << endl;
        return 1;
    }

    ifstream parse_in(argv[1]);
    if (!parse_in.is_open())
    {
        cerr << "Failed to open input file: " << argv[1] << endl;
        return 1;
    }
    vector<logic_info> logic_blocks_table;
    vector<pin_info> pin_table;
    vector<net_info> net_table;
    unordered_map<string, int> name2idx;
    unordered_map<string, int> pin2idx;
    unordered_map<string, int> net2idx;
    unordered_map<string, net_hpwl_cong> net2hpwl_cong;
    double total_hpwl, CC;
    parse_input(parse_in, logic_blocks_table, pin_table, net_table, name2idx, pin2idx, net2idx);
    vector<vector<string>> array2D(R, vector<string>(C, "none"));
    vector<vector<int>> U(R, vector<int>(C, 0));

    initial_place(array2D, R, C, logic_blocks_table);

    compute_hpwl_cong(logic_blocks_table, pin_table, net_table, U, net2hpwl_cong, total_hpwl, CC);
    SA(logic_blocks_table, array2D,
       pin_table, net_table,
       name2idx, pin2idx,
       net2hpwl_cong, U,
       total_hpwl, CC);

    cout << "[FINAL] HPWL=" << total_hpwl
         << " CC=" << CC << endl;
    write_output(string(argv[2]), logic_blocks_table);

    auto end = Clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    cout << "Execution time: " << duration.count() << " milliseconds" << endl;

    return 0;
}
