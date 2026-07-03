function fs_teardown()
%FS_TEARDOWN  Standard test teardown: print summary, exit with 0 or 1.
%
%   Reads appdata slot 'fs_test_state'; prints pass/fail counts + timing;
%   exits process with code 0 (all pass) or 1 (any fail).

    st = getappdata(0, 'fs_test_state');
    if isempty(st), st = struct('pass', 0, 'fail', 0, 'first_fail', ''); end

    name = getappdata(0, 'fs_test_name');
    if isempty(name), name = '(unnamed)'; end

    t0 = getappdata(0, 'fs_test_t0');
    if isempty(t0), elapsed = 0; else, elapsed = toc(t0); end

    total = st.pass + st.fail;
    fprintf('\n──────────────────────────────────────────────────\n');
    if st.fail == 0 && total > 0
        fprintf('  TEST OK    %s   %d/%d passed   (%.2fs)\n', name, st.pass, total, elapsed);
        exit(0);
    elseif total == 0
        fprintf('  TEST WARN  %s   no assertions recorded   (%.2fs)\n', name, elapsed);
        exit(0);
    else
        fprintf('  TEST FAIL  %s   %d/%d passed   %d failed   (%.2fs)\n', ...
                name, st.pass, total, st.fail, elapsed);
        fprintf('  first failure: %s\n', st.first_fail);
        exit(1);
    end
end
